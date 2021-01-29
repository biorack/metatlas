with open('/global/homes/b/bpb/repos/metatlas/metatlas/tools/20200306_KBL-AK_GM_504788_Blumenol_Final_QE-HF_C18-CustMog_USDAY48854_worklist.sh','r') as fid:
    j = fid.read()
counter = 0
nj = []
for i,jj in enumerate(j.split('\n')):
    nj.append('%s &'%jj)
    counter += 1
    if counter == 10:
        nj.append('wait')
        counter = 0
with open('/global/homes/b/bpb/repos/metatlas/metatlas/tools/20200306_KBL-AK_GM_504788_Blumenol_Final_QE-HF_C18-CustMog_USDAY48854_worklist_bundles.sh','w') as fid:
    fid.write('\n'.join(nj))
