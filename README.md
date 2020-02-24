# CADPlaneFitting
Fit Planes for CAD models

![Plane Fitting Results](https://github.com/hjwdzh/CADPlaneFitting/raw/master/res/teaser.png)

### Compile
```
mkdir build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

### Run
Take a CAD model as input. Output a list of Plane IDs for each face. Optionally output a colored model with different planes.
```
./plane input.obj output.txt [visual.obj]
```