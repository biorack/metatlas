
# Simple workflow to make your own atlas

```
my_mz_ref = metob.MzReference(mz=1000)
my_rt_ref = metob.RtReference(rt_peak=3)
my_compound = metob.Compound()
my_id = metob.CompoundIdentification(compound=[my_compound],rt_reference=[my_rt_ref],mz_reference=[my_mz_ref])
my_atlas = metob.Atlas(name='my first atlas',compound_identifications=[my_id])
```
