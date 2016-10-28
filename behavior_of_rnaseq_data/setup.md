## Setup

There are now a number of different components to this so it's worth describing my set up. In case it's not clear, all of these steps are meant to be followed in the current directory.

## Creating a dedicated and isolated environment
[Install miniconda](http://conda.pydata.org/miniconda.html) if you haven't already done so.

```shell
# Create "immuneinf" python3 environment
$ conda create -n immuneinf python=3 matplotlib pandas seaborn numpy scipy jupyter
```

## Activating the environment and updating the dependencies
Once you have your conda environment ready,
to activate it, use:

```shell
$ source activate immuneinf
```

and to update the dependencies, use:
```shell
(immuneinf) $ pip install -r requirements.txt
```

To add this kernel to your jupyter environment (e.g. to be able to select it for notebooks), run:

```shell
(immuneinf) $ python -m ipykernel install --user --name=immune3
```

This will add a new kernel named "immune3" to your jupyter environment, referencing the "immuneinf" conda environment.

## re-execute the notebooks

Finally, given the amount of time it takes for these models to fit (currently!), I typically run them overnight.

To set this up, first install `nbutils` from github:

```shell
(immuneinf) $ pip install git+git://github.com/jburos/nbutils
```

This will add an executable to your environment, called `nbexecute`.

You can then re-execute a notebook (or several notebooks) as follows:

```shell
(immuneinf) $ nbexecute --kernel-name=immune3 --timeout 900000 0.85*.ipynb
```

The line above will execute each of the `0.85*.ipynb` files in sequential order, referencing the `immune3` kernel (this technically isn't required but I find it helpful to be explicit -- by default `nbexecute` will use the kernel from the last checkpoint).

What's cool about this, is that **the models are cached** (hooray!). This means that, the next day, you can re-execute the notebook and it will be reasonably fast.

