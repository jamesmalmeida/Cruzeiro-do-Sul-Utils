# Cruzeiro-do-Sul-Utils
Libraries for Cruzeiro do Sul Database


# XDI -> json example:

```
import CZDS_utils
from CZDS_utils import XDI

#instantiate XDI class
xdi = XDI("somefile.xdi")

# print json string before normalization
print(xdi)

# include normalized data column 
xdi.normalize()

#print json string after normalization
print(xdi)

#save to file (json)

xdi.save("somefile.json")
```
