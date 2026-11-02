#!/usr/bin/env python3
"""Convert copied upstream numerical input inventory to independently testable CSV."""
import csv, json
from pathlib import Path
root=Path(__file__).resolve().parent
with (root.parent/'fixtures/catalog_inputs.csv').open('w',newline='') as f:
 writer=csv.writer(f)
 writer.writerow(['id','func','a','b','c','z'])
 for i,line in enumerate((root/'sci_special_cases_hypergeometric.jsonl').read_text().splitlines()):
  case=json.loads(line)
  name,args=case['name'],case['args']
  if name=='hyp0f1': a,b,c,z='',args[0],'',args[1]
  elif name=='hyp2f1': a,b,c,z=args
  elif name=='wrightomega': a,b,c,z='','','',args[0]
  else:a,b,c,z=args[0],args[1],'',args[2]
  writer.writerow([f'S{i:04d}',name,a,b,c,z])
