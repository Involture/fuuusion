from ops import *
import scipy.misc

IMAGENAME = "roller"

binm = np.load("res/" + IMAGENAME + "_bin-.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_bin-.jpg", binm)
bin = np.load("res/" + IMAGENAME + "_bin.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_bin.jpg", bin)
cbin = np.load("res/" + IMAGENAME + "_closedbin.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_closedbin.jpg", cbin)
contour = np.load("res/" + IMAGENAME + "_contour.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_contour.jpg", contour)
fft = sf(np.abs(np.load("res/" + IMAGENAME + "_fft.npy")), 1000)
scipy.misc.imsave("res/" + IMAGENAME + "_fft.jpg", fft)
filtered = np.load("res/" + IMAGENAME + "_filtered.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_filtered.jpg", filtered)
lab = np.load("res/" + IMAGENAME + "_lab.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_lab.jpg", lab)
localised = np.load("res/" + IMAGENAME + "_localised.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_localised.jpg", localised)
resp = np.load("res/" + IMAGENAME + "_resp.npy")
scipy.misc.imsave("res/" + IMAGENAME + "_resp.jpg", resp)
