#encoding=utf-8  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
 
#自定义函数 e指数形式
def func(x,lb,c,d0,const):
    return -c*(1+(x-d0)/lb)*np.exp(-(x-d0)/lb)+const
 
x = []
y = []
#定义x、y散点坐标
with open('wat') as f:
    for line in f:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

x=np.array(x)
y=np.array(y)
 
#非线性最小二乘法拟合
popt, pcov = curve_fit(func, x, y)
#获取popt里面是拟合系数
lb = popt[0]
c = popt[1]
d0 = popt[2]
const=popt[3]
yvals = func(x,lb,c,d0,const) #拟合y值
print(lb,c,d0)
 
#绘图
xrange=np.linspace(1.2,10.3,100)
plot1 = plt.plot(x, y, 's',label='original values')
plot2 = plt.plot(xrange, func(xrange,lb,c,d0,const), 'r',label='polyfit values')
plt.xlabel('x')
plt.ylabel('y')
plt.title('LLZO-Li UBER')
plt.show()
plt.savefig('test2.png')

