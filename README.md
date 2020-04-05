### 需要手动调整的参数

* 15-17行、19-21行、23-25行处的 n, bestRoute, dir 

不同的实例请采用我的设置，15-17行对应Sahara，19-21行对应Djibouti，23-25行对应Qatar

n 为城市数目，bestRoute为全局最优值（注 national TSP计算城市之间的最短距离时取四舍五入的整数，导致了与我设置的全局最优值有误差），dir为输入文件的相对路径

输入数据格式需要与此仓库下的.txt文件保持一致或者直接采用此仓库下的数据

如果顺利结果应该和图中类似

![](./qatar_result.png)
![](./sahara_result.png)
![](./djibouti_result.png)

### 可选调整参数

* 27-33行处的参数，可参考注释
