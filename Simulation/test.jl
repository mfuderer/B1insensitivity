using PythonPlot

_,ax = PythonPlot.subplots(1); fig = gcf()
bbb = rand(1000)
ax.plot(bbb)
fig.suptitle("bbb")
PythonPlot.show()

aaa = rand(100)
figure()
plot(aaa) 
title("aaa")
PythonPlot.show()
