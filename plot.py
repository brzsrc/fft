import matplotlib.pyplot as plt
x_value = [2,4,8,16,32,64,128,256]
# y_value = [0.4505,0.7022,0.9465,1.1133,1.1128,1.1128,1.1128,1.1128]
y_value = [1225425882.3295,
1067724913.2767,
988715429.1227,
948331084.0359,
957531963.2041,
976005209.4008,
1013280002.3975,
1082031806.6812]
plt.figure(figsize=(10,7))
plt.plot(x_value,y_value)
plt.scatter(x_value,y_value)
plt.ylabel('the total energy required',size=15)
# plt.ylabel('sim_IPC',size=15)
plt.xlabel('LSQ size',size=15)
# plt.title('？',size=20)
plt.show()