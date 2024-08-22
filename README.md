# Cruzeiro-do-Sul-Utils
Libraries for Cruzeiro do Sul Database


# XDI -> json example:

```
import CZDS_utils
from CZDS_utils import XDI

xdi = XDI("somefile.xdi")
print(xdi) #print json string before normalization
xdi.normalize() #include normalized data column 
print(xdi) #print json string after normalization
```
