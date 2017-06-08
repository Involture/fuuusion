import time


n = 200000
t0 = time.time()
for i in range(n):
	print(i)
print((time.time() - t0) / n)

