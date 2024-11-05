# transnetyx-pyrat-utils

Export your transnetyx genotyping results as a .csv file.

## Python

Just supply the name of the file you have in your Downloads folder, and run ``.convert()``

```python
from transnetyx_pyrat_utils import Converter
converter = Converter("OrderResults-2215183-638663954144327998.csv").convert()
```

## CLI 

```bash
pdm run transpyrat OrderResults-2215183-638663954144327998.csv
```

or 

```bash
pdm run transpyrat -f OrderResults-2215183-638663954144327998.csv
```

or 

```bash
pdm run transpyrat --file OrderResults-2215183-638663954144327998.csv
```

get the full help with

```bash
pdm run transpyrat --help
```
⤵️
```bash
Usage: -c [-h] [-f FILE] [file]

Positional Arguments:
  file                  The file name or path leading to the file to be converted. By default, the file name is searched for    
                        in the Downloads folder of the current user's home.

Options:
  -h, --help            show this help message and exit
  -f, --file FILE       The file name or path leading to the file to be converted. By default, the file name is searched for    
                        in the Downloads folder of the current user's home.
```

## Config

As you have more lines and strains to convert from transnetyx framework intopyrat one, add these strains in the `transnetyx_config.json` configuration file, that sits in the source code folder of this package. (it can also be supplied from a different path, externally, by changing the transnetyx_config_filename or pyrat_config_filename attributes of your converter object, after instanciating it with a filename, and before running ``.convert()``)