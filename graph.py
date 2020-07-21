# импортируем модули
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
# уравнение поверхности
f = lambda x, y: x ** 2 - y ** 2
# создаём полотно для рисунка
fig = plt.figure(figsize = (8, 6))
# создаём рисунок пространства с поверхностью
ax = fig.add_subplot(1, 1, 1, projection = '3d')
# размечаем границы осей для аргументов
xval = np.linspace(-4, 4, 100)
yval = np.linspace(-4, 4, 100)
# создаём массив с xval столбцами и yval строками
# - в этом массиве будут храниться значения z
x, y = np.meshgrid(xval, yval)
# приравниваем z к функции от x и y 
z = f(x, y)
print(z)
# создаём поверхность
surf = ax.plot_surface(
# отмечаем аргументы и уравнение поверхности
x, y, z, 
# шаг прорисовки сетки
# - чем меньше значение, тем плавнее
# - будет градиент на поверхности
rstride = 2,
cstride = 2,
# цветовая схема plasma
cmap = cm.viridis)
plt.show()

#k=input("press close to exit") 