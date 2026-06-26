title: Monitors
---

[TOC]

The monitors are specified individually as blocks in the control file.
The general keywords that can be specified are explained in the table below:

| Keyword                 | Description                                                                                   | Default value |
|-------------------------|-----------------------------------------------------------------------------------------------|---------------|
| monitors flush interval | *INTEGER*: Iteration interval to flush the monitor information to the monitor files.          | 100           |


## Residual Monitors

## Statistics Monitor 

```markdown
#define statistics
   initial time      = 1.d0
   initial iteration = 10
   sampling interval = 10
   dump interval     = 20
   @start
#end
```

By default, the statistic monitor will average following variables:

<div class="multicol">
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>u</li>
            <li>v</li>
            <li>w</li>
        </ul>
    </div>
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>uu</li>
            <li>vv</li>
            <li>ww</li>
        </ul>
    </div>
    <div style="float: left; width: 33.33%;">
        <ul>
            <li>uv</li>
            <li>uw</li>
            <li>vw</li>
        </ul>
    </div>
    <div style="clear: both;"></div>
</div>



A keyword preceded by @ is used in real-time to signalize the solver what it must do with the statistics computation:

<div class="multicol">
    <ul>
        <li>@start</li>
        <li>@pause</li>
        <li>@stop</li>
        <li>@reset</li>
        <li>@dump</li>
    </ul>
</div>

After reading the keyword, the solver performs the desired action and marks it with a star, e.g. @start*.

@note 
Real-time keywords may not work in parallel MPI computations. It depends on how the system is configured.
@endnote

## Probes

```markdown
#define probe 1
   name     = SomeName
   variable = SomeVariable
   position = [0.d0, 0.d0, 0.d0]
#end
```

| Keyword   | Description                                                 | Default value       |
|-----------|-------------------------------------------------------------|---------------------|
| name      | *CHARACTER*: Name of the monitor.                           | **Mandatory Keyword** |
| variable  | *CHARACTER*: Variable to be monitored. Implemented options are: pressure, velocity, u, v, w, mach, k. | **Mandatory Keyword** |
| position  | *REAL(3)*: Coordinates of the point to be monitored.        | **Mandatory Keyword** |

### Probes from a file

Instead of (or in addition to) defining probes one-by-one with `#define probe` blocks,
a set of probes can be defined in bulk through a plain text `.dat` file, each line
holding the coordinates of a point. The same list of variables is sampled at every
probe in the file, and the probes-file block is defined like any other monitor block:

```markdown
#define probe file
   file      = MyProbes.dat
   variables = pressure u v
#end
```

| Keyword   | Description                                                          | Default value          |
|-----------|-----------------------------------------------------------------------|-------------------------|
| file      | *CHARACTER*: Path to the `.dat` file defining the probes' coordinates. | **Mandatory Keyword** |
| variables | *CHARACTER*: Whitespace-separated list of variables to sample at every probe in the file. | **Mandatory Keyword** |

```
# x       y      z
0.5d0     0.d0   0.d0
1.0d0     0.d0   0.d0
2.0d0     1.d0   0.d0
```

Lines that are empty or start with `#` are ignored; every remaining line holds the
`x y z` coordinates of one probe. For each of these probes a dedicated output file
`<solution_file>.probe_N.probe` is created (`N` being the probe's position in the
file), with one column per variable (as listed in `variables`) and one row per saved
time step (in addition to the `Iteration` and `Time` columns).

To keep the per-iteration screen log readable when many probes are defined this way,
file-based probes are excluded from it; their individual `.probe` files are still
written normally. Instead, a summary of the probes-file monitor (file path, number of
probes and sampled variables) is printed once in the startup log, right after the
"Time integrator" section.


## Surface Monitors

```markdown
#define surface monitor 1
   name              = SomeName
   marker            = NameOfBoundary
   variable          = SomeVariable
   reference surface = 1.d0
   direction         = [1.d0, 0.d0, 0.d0]
#end
```
| Keyword           | Description                                                                                      | Default value       |
|-------------------|--------------------------------------------------------------------------------------------------|---------------------|
| name              | *CHARACTER*: Name of the monitor.                                                                | **Mandatory Keyword** |
| marker            | *CHARACTER*: Name of the boundary where a variable will be monitored.                             | **Mandatory Keyword** |
| variable          | *CHARACTER*: Variable to be monitored. Implemented options are: mass-flow, flow, pressure-force, viscous-force, force, lift, drag, pressure-average. | **Mandatory Keyword** |
| reference surface | *REAL*: Reference surface [area] for the monitor. Needed for "lift" and "drag" computations.    | --                  |
| direction         | *REAL(3)*: Direction in which the force is going to be measured. Needed for "pressure-force", "viscous-force" and "force". Can be specified for "lift" (default [0.d0,1.d0,0.d0]) and "drag" (default [1.d0,0.d0,0.d0]). | --                  |

## Volume monitors

Volume monitors compute the average of a quantity in the whole domain. They can be scalars(s) or vectors(v).

```markdown
#define volume monitor 1
   name     = SomeName
   variable = SomeVariable
#end
```

| Keyword  | Description                                                                                                    | Default value       |
|----------|----------------------------------------------------------------------------------------------------------------|---------------------|
| name     | *CHARACTER*: Name of the monitor.                                                                              | **Mandatory Keyword** |
| variable | *CHARACTER*: Variable to be monitored. The variable can be scalar (s) or vectorial (v). Implemented options are: kinetic energy (s), kinetic energy rate (s), enstrophy (s), entropy (s), entropy rate (s), mean velocity (s), velocity (v), momentum (v), source (v). | **Mandatory Keyword** |


## Load Balancing Monitors

Load balancing monitors compute the DOF in each partition of the mesh to check the unbalance. 

```markdown
#define load balancing monitor 1
   name     = SomeName
   variable = SomeVariable
#end
```

| Keyword  | Description                                                                                                     | Default value       |
|----------|-----------------------------------------------------------------------------------------------------------------|---------------------|
| name     | *CHARACTER*: Name of the monitor.                                                                               | **Mandatory Keyword** |
| variable | *CHARACTER*: Variable to be monitored. Implemented options are: max dof per partition, min dof per partition, avg dof per partition, absolute dof unbalancing, relative dof unbalancing. | **Mandatory Keyword** |


The variable *absolute dof unbalancing* is computed as follows:

\[
\text{absolute dof unbalancing} = \text{max dof per partition} - \text{min dof per partition}
\]

The variable *relative dof unbalancing* is computed as follows:

\[
\text{relative dof unbalancing} = \frac{100 \times \text{absolute dof unbalancing}}{\text{avg dof per partition}}
\]


