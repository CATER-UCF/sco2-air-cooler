# Heat Exchanger 2D

The goal of this project is to create a two-dimensional, transient model of a countercurrent crossflow heat exchanger. Here, we use a flow network of zero-dimensional heat exchanger elements.

Currently, there are two different network models. With mixers and separators, the air streams are combined at the inlet of each tube pass. Without mixing, there is a separate air stream for each column of elements. Show below are networks with three passes and four elements per pass. The actual flowsheet model has eight and seven respectively.

### Heat Exchanger Network with Mixers and Separators

![](./images/2D_network_with_mixing.PNG)

### Heat Exchanger Network without Mixing

![](./images/2D_network_no_mixing.PNG)


## Code Gen

To the best of my knowledge, we cannot have a list or array of unit models added to an IDAES flowsheet. So each heat exchanger element needs to be uniquely named. To handle this, I wrote a program that generates code blocks that need to be pasted into a flowsheet:

```
python ./flowsheets/code_gen.py
```

Results are added to [here](./flowsheets/sco2_2d_steady_state). The code gen program can work for any number of heat exchanger passes and elements per pass. By default, heat exchanger elements are named `e0`, `e1`, `e2`, ... and so forth.

## Status / Issues

The 2D model is showing reasonable results for steady-state. Without air stream mixing, the flowsheet solves with little difficulty. With mixing, it becomes harder to solve. In any case, the next steps are to:

1. Add the constraints and equations to make HTC's a function of fluid properties. 
2. Add a lumped capacitance term to the 0D model's energy balance. This is started [here](./models/heat_exchanger_0d_dynamic.py).

However... I'm having difficulty setting up a dynamic flowheet. Here, I tried using the standard 0D IDAES heat exchanger as a baseline:

```
python ./flowsheets/sco2_0d_dynamic.py
```

The model solves but the results don't make sense yet. Inlet conditions are fixed at every time step. The outlet mass flow flows, though, don't match the inlet mass flows. I tried fixing the outlet mass flows but then the model doesn't solve.
