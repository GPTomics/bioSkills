'''Load an FCS file in Python: FlowKit for analysis, readfcs for the scanpy/AnnData bridge.'''
# Reference: flowkit 1.1+, readfcs 1.1+ | Verify API if version differs
import glob
import flowkit as fk
import readfcs

fcs_files = sorted(glob.glob('data/*.fcs'))
print(f'found {len(fcs_files)} FCS files')

# FlowKit: full event/metadata access; source can be raw, comp, or xform
sample = fk.Sample(fcs_files[0])
print('channels:', sample.pnn_labels)        # detector names
print('antibodies:', sample.pns_labels)      # $PnS antibody descriptions
events = sample.as_dataframe(source='raw')
print('event matrix:', events.shape)

# readfcs: FCS -> AnnData for the scanpy/scverse ecosystem
adata = readfcs.read(fcs_files[0])
print('AnnData:', adata.shape)                # cells x channels
print(adata.var.head())                       # channel + antibody mapping
