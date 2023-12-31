import numpy as np
import cv2
import csv

write_result = True

N = 500
template = np.zeros((N,N), dtype = np.uint8)
texts = ["Happy", "New", "Year"]
font = cv2.FONT_HERSHEY_SIMPLEX
fontscale = 4.0
color = 255
thickness = 15

sizes = [cv2.getTextSize(text, font, fontscale, thickness)[0] for text in texts]
y = 150

for text, size in zip(texts, sizes):
    x = (template.shape[1] - size[0])//2
    org = (x, y)
    print(org)
    template = cv2.putText(template, text, org, font, fontscale, color, thickness, cv2.LINE_AA, False)
    y += size[1] + 50

where_text = (template != 0)
canvas = np.zeros((N,N), dtype = np.uint8)
n_disks_total = 400
n_disks = 0
radius_min = 3
radius_max = 8
centers = []
radii = []

while n_disks < n_disks_total:
    where_not_disks = (canvas == 0)
    candidates_y, candidates_x = np.where(where_text & where_not_disks)
    index = np.random.randint(len(candidates_y))
    center = np.array([candidates_x[index], candidates_y[index]])

    if n_disks < 1:
        radius = radius_max
    elif n_disks == n_disks_total:
        break
    else:
        radius = min(np.linalg.norm(center - center_other) - radius_other for center_other, radius_other in zip(centers, radii))
        radius = min(radius, radius_max)
        radius = int(radius)

    if radius >= radius_min:
        centers.append(center)
        radii.append(radius)
        n_disks += 1

        canvas = cv2.circle(canvas, center, radius, color, -1)
        print(f"{n_disks}/{n_disks_total}")
    

cv2.imshow("", canvas)
cv2.waitKey(0)
cv2.destroyAllWindows()

if write_result:
    with open(r"C:\Users\bart1\Documents\Julia_projects\Rays\experiments\happy_new_year\disks.csv", 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["center_x", "center_y", "radius"])
        for i in range(n_disks_total):
            center = centers[i]
            radius = radii[i]
            writer.writerow([*(center/N)[::-1], radius/N])
