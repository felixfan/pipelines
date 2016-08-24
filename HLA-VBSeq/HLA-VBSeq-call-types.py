import sys

def call_hla_types(iid, doc):
      '''
      call HLA types from HLA-VBSeq outputs     
      iid: individual id, there are 6 files for each individual    
      doc: threshold of depth of coverage, 20% of the depth of coverage was recommended     
      '''
      ans = [iid]
      genes = ['HLA_A','HLA_B','HLA_C','HLA_DQA1','HLA_DQB1','HLA_DRB1']
      for gene in genes:
            fn = iid + '.' + gene + '.txt'
            f = open(fn)
            n = 0
            for line in f:
                  n += 1
                  if n == 1:
                        line = line.strip()
                        ws = line.split()
                        print ws[1]
                        if float(ws[1]) < doc:
                              ans.append('NA')
                              ans.append('NA')
                              break
                        else:
                              g1 = ws[0]
                              c1 = float(ws[1])
                  elif n == 2:
                        line = line.strip()
                        ws = line.split()
                        if float(ws[1]) < doc:
                              if c1 >= 2.0 * doc:
                                    ans.append(g1)
                                    ans.append(g1)
                              else:
                                    # heter?????????????
                                    ans.append(g1)
                                    ans.append(ws[0])
                        else:
                              if c1 >= 2.0 * float(ws[1]):
                                    ans.append(g1)
                                    ans.append(g1)
                              else:
                                    ans.append(g1)
                                    ans.append(ws[0])
                  else:
                        break
            f.close()
      for a in ans:
            print a,
      print

if __name__ == '__main__':
      iid = sys.argv[1]
      doc = float(sys.argv[2])
      call_hla_types(iid, doc)
