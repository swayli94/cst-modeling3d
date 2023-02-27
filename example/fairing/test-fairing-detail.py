

import numpy as np
from cst_modeling.surface import BasicSurface

from scipy.interpolate import CubicSpline, CubicHermiteSpline
import matplotlib.pyplot as plt

#* 插值方法
# https://docs.scipy.org/doc/scipy/reference/interpolate.html


if __name__ == "__main__":

    #*==============================================
    # Step 1: 确定鼓包的总体参数 in Fairing.txt
    #   1. 鼓包是由若干 垂直于流向的截面曲线 拉伸生成的
    #   2. 总体参数：   各个截面的流向位置（控制文件每一行中的第三个参数）
    #!  需要定义的量：n_sec, Fairing.txt 每一行中的第三个参数, nn (等于tube-curve.dat中的点数)
    #*==============================================

    fairing = BasicSurface(n_sec=8, name='Fairing', nn=51, ns=101, projection=False)
    fairing.read_setting('Fairing.txt')

    #*==============================================
    # Step 2: 读入鼓包与机身的过渡位置的流向截面曲线 tube-curve.dat
    #   1. 从 CATIA 中截取该曲线，整理为 Tecplot 格式
    #   2. 读入该曲线中鼓包所对应范围的曲线 （本例子中，该曲线为相位角 [90deg, 180deg] 的一段）
    #      则该曲线为鼓包各个流向控制截面的基础型线，
    #      鼓包的控制参数定义了在基础型线上的几何叠加量
    #   3. 将之转换为单位曲线后，赋值给鼓包的各个控制曲面
    #      注：转换为单位曲线是因为程序的定义
    #!  需要定义的量：tube-curve.dat
    #*==============================================

    # Read raw data
    with open('tube-curve.dat', 'r') as f:
        X_end = []
        Y_end = []
        lines = f.readlines()
        for i in range(len(lines)-2):
            line = lines[i+2].split()
            X_end.append(float(line[0]))
            Y_end.append(float(line[1]))

    X_end = np.array(X_end)
    Y_end = np.array(Y_end)

    for i in range(fairing.n_sec):
        fairing.secs[i].xx = X_end.copy()
        fairing.secs[i].yy = Y_end.copy()

    #*==============================================
    # Step 3: 定义鼓包控制截面的几何叠加量
    #   1. 鼓包控制截面的最终几何 = 机身的流向截面曲线 + 三次样条定义的几何叠加量
    #   2. 三次样条由若干控制点来定义，使用绝对坐标
    #!  需要定义的量：各个控制截面的几何叠加量的控制参数
    #*==============================================

    # 三次样条的控制点
    X0 = [X_end[-1], -2.20, -1.50, -0.80, X_end[0]]
    Y0 = [0.0,        0.80,  0.60,  0.10, 0.0]
    curve = CubicSpline(X0, Y0, bc_type=((1,0.0), (1,0.0)))

    '''
    X0 =   [X_end[-1], -2.20, -1.20, X_end[0]]
    Y0 =   [0.0,        0.80,  0.30, 0.0]
    dydx = [0.0,        0.0,  -0.2,  0.0]
    curve = CubicHermiteSpline(X0, Y0, dydx=dydx)
    '''

    y_increment = curve(X_end)

    if False:
        #! Plot curves on screen
        plt.plot(X_end, Y_end, 'r')
        plt.plot(X0, Y0, 'go')
        plt.plot(X_end, y_increment, 'g')
        plt.plot(X_end, Y_end+y_increment, 'b')
        plt.axis('equal')
        plt.legend(['baseline curve', 'control points', 'incremental curve', 'final curve'])
        plt.show()

    # 将几何增量 y_increment 叠加到原始曲线上
    fairing.secs[2].yy += y_increment
    fairing.secs[3].yy += y_increment
    fairing.secs[4].yy += y_increment*1.2
    fairing.secs[5].yy += y_increment*1.2

    #*==============================================
    # Step 4: 生成三维曲面，光滑曲面，输出
    #!  需要定义的量：smooth 光滑范围的截面编号
    #*==============================================
    fairing.geo()
    
    fairing.smooth(1, 3, smooth0=True, smooth1=True)
    fairing.smooth(5, 7, smooth0=True, smooth1=True)

    fairing.flip(axis='+Z +Y')

    fairing.output_tecplot(fname='fairing.dat')

    fairing.output_plot3d(fname='fairing.grd')

