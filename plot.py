import matplotlib.pyplot as plt


def draw_plot(x, y):
    plt.plot(x, y)
    plt.title("Electric potential equation")
    plt.xlabel("x")
    plt.ylabel("FEM result")
    plt.show()
