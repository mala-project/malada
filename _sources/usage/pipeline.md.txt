# Pipeline

The central class of MALADA workflows is the `DataPipeline` class. It (semi-automates) the entire MALADA workflow (see [current limitations](global_parameters.md))
To use it, create a `Parameters` object, assign the parameters you want, create a `Pipeline` object with it, and run:

```python
# Create some parameters and fill them with life.
params = malada.Parameters()
params.number_of_atoms = 2
params.temperature = 298
params.element = "Be"

# Create a pipeline and let it run.
pipeline = malada.DataPipeline(params)
pipeline.run()
```

## Providers

Internally, the MALADA pipeline works with "Providers". Each such a class has a specific task, and specific in/outputs. Each provider uses the output of the preceding step to generate input for the next one. The pipeline can be customized by giving it a customized provider upon creation:


```python
# Create some parameters and fill them with life.
params = malada.Parameters()
params.number_of_atoms = 2
params.temperature = 298
params.element = "Be"

# Create a custom CrystalStructureProvider, because we already have a crystal
# structure downloaded to file.
crystal_structure = malada.CrystalStructureProvider(params,
                                                    external_cif_file=
                                                    "path/to/some/file")

# Create a pipeline and let it run.
pipeline = malada.DataPipeline(params, crystal_structure_provider=crystal_structure)
pipeline.run()
```

This feature is especially helpful given that currently, some automation features are not yet implemented in MALADA.
The providers can also be used in standalone to perform only a specific processing step. For this, such a provider object has to be called, and then executed via the `provide()` function.

```python 
params = malada.Parameters()
params.number_of_atoms = 2
params.temperature = 298
params.element = "Be"

# Run a custom MD calculation.
md = malada.MDProvider(params)
md.provide("path/to/runfolder", "some_supercell",
           "some_convergence_parameters.xml",
            "some_parallelization_info.md")
```

This piece of code will execute and analyse an MD calculation in the specified folder, using a specific supercell, DFT parameters and parallelization strategy. 
Which inputs have to be given to the `provide` function differs. For a full list of providers, see the [list of providers](list_of_providers.md).

## Runners 

MALADA does NOT aim to interface with DFT/MD codes directly. It is NOT supposed to perform tasks similar to e.g. ASE, which can directly run DFT/MD calculations. Rather it just provides inputs, that can then be run with some run system (e.g. via slurm or ASE). See [current limitations](global_parameters.md) for an overview of what is currently implemented in that regard.
