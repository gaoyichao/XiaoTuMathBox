import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('bienao.csv')

plt.figure(figsize=(10, 6))
plt.style.use('seaborn-v0_8-whitegrid')


# 4. 根据要求分别绘制不同的曲线
x = df['x']

# y: 实线，实心圆标记
plt.plot(x, df['y'], linestyle='-', marker='o', color='black', label='y', linewidth=2)

# p5, p7, p9: 虚线，无标记
plt.plot(x, df['p5'], linestyle='--', marker='None', label='p5', linewidth=2)
plt.plot(x, df['p7'], linestyle='--', marker='None', label='p7', linewidth=2)
plt.plot(x, df['p9'], linestyle='--', marker='None', label='p9', linewidth=2)

# p15: 虚线，空心方框标记
plt.plot(x, df['p15'], linestyle='--', marker='s', markerfacecolor='none', markeredgecolor='auto', label='p15', linewidth=2)

# p17: 虚线，空心三角标记
plt.plot(x, df['p17'], linestyle='--', marker='^', markerfacecolor='none', markeredgecolor='auto', label='p17', linewidth=2)


# 5. 添加图表细节
plt.title('CSV Data Curves', fontsize=16)
plt.xlabel('X Axis', fontsize=14)
plt.ylabel('Values', fontsize=14)
plt.legend(title='Parameters', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)

# 6. 显示图表
plt.tight_layout()
plt.show()
