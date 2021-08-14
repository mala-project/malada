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
The providers can also be used in standalone to perform only a specific processing step 

```python 
params = malada.Parameters()
params.number_of_atoms = 2
params.temperature = 298
params.element = "Be"

md = malada.MDProvider(params)
md.provide("path/to/folder", "some_supercell",
           "some_convergence_parameters.xml",
            "some_parallelization_info.md")
```

This piece of code will execute only an MD calculation in the specified folder. 
For a full list of providers, see the [list of providers](list_of_providers.md).
