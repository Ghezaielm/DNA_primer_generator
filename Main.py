class DNA_primer_generator(): 
    
    def __init__(self):
        
        self._sequences = False
        self._features = False
        self.sequence = False
        self.features = False
        self.DNA_base = {"A":1,"T":2,"G":3,"C":4}
        self.RNA_base = {"A":1,"U":2,"G":3,"C":4}
        self.DNA_comp = {"A":"T","T":"A","G":"C","C":"G"}
        self.RNA_comp = {"A":"U","U":"A","G":"C","C":"G"}
    
    # Returns reverse complementary
    def reverseComp(self,s,current="Y"):
        if current in ["Y","y","Yes","yes","YES"]:
            if "U" not in s:
                s = [self.DNA_comp[i.upper()] for i in s]
                return s
            else:
                s = [self.RNA_comp[i.upper()] for i in s]
                return s 
        else:
            return False
        
    # To handle N containing sequences
    def setNomenclature(self,dic,_type="DNA"):
        if _type=="DNA":
            for i in dic:
                self.DNA_base[i]=dic[i]
                self.DNA_comp[i]=dic[i]
        elif len(dic)<1:
            return False
        else:
            self.RNA_base[i]=dic[i]
            self.RNA_comp[i]=dic[i]
    
    # Fetch nucleotide frequencies
    def nucleotideFrequencies(self,s):
        return {j:sum([1 if i.upper()==j else 0 for i in s ]) for j in s}
    
    # Get sequence length
    def getLength(self):
        return len(self.current)
    
    # Get melting temperature
    def getTm(self,s):
        freqs = DNA_analyzer.nucleotideFrequencies(self,s)
        for j in ["A","T","G","C"]:
            if j not in freqs:
                freqs[j]=0
        if len(s)<13:
            return 4*(freqs["G"]+freqs["C"])+2*(freqs["A"]+freqs["T"])
        else:
            return 64.9 + 41*(freqs["G"]+freqs["C"]-16.4)/len(s)
        
    # Measure complementarity between sequences
    def measureComp(self,a,b):
        affinity = 0
        score = 0
        b = DNA_analyzer.reverseComp(self,b)
        for i in range(len(a)):
            if a[i]==b[i]:
                score+=1
            if i!=len(a)-1: 
                if a[i+1]==b[i+1]:
                    affinity+=1

            
        if affinity == 0 and score>0:
            return score
        else:
            return affinity*score
        
    # Find primers on a set of sequences
    def FindPrimers(self,s,start,end,max_size=26,margins=10,tmin=52,tmax=58):
        
        ### Check if a primer is a palindrom
        def isPalindrom(s):
            s = [self.DNA_base[i.upper()] for i in s]
            
            def gamma(a,b):
                return((a+2)/a)*b
            if functools.reduce(gamma,s)==functools.reduce(gamma,s[::-1]):
                return True
            else: 
                return False   
        
        # Right and left sides of the primers
        if len(s)<max_size:
            max_size = round(len(s)/2)
        max_size = round(max_size/2)
        
        # Convert sequences to digits
        strand = s
        
        # Placeholders for putative primers
        direct_primers = []
        indirect_primers = []
        
        # Scan the direct strand
        for idx in range(3,max_size):
            for left in range(start-margins,start):
                primer = DNA_analyzer.reverseComp(self,strand[left-idx:left+idx])
                if len(primer)==0:
                    pass
                elif isPalindrom(primer):
                    pass
                else:
                    tm = DNA_analyzer.getTm(self,primer)
                    if tm>=tmin and tm<=tmax:
                        direct_primers.append((left,idx,tm))
                        
        # Scan the reverse strand                
        for idx in range(3,max_size):
            for right in range(end,end+margins):
                primer = strand[right-idx:right+idx]
                if len(primer)==0:
                    pass
                elif isPalindrom(primer):
                    pass
                else:
                    tm = DNA_analyzer.getTm(self,primer)
                    if tm>=tmin and tm<=tmax:
                        indirect_primers.append((right,idx,tm))
        self.direct_primers = [i for i in direct_primers if i[0]>1]
        self.indirect_primers = [i for i in indirect_primers if i[0]>1]
        self.strand = strand
        
    # Check primers specificity     
    def CheckPrimers(self):
        # Test primers for specificity:
        direct_primers = self.direct_primers
        indirect_primers = self.indirect_primers
        strand = self.strand
        direct_primers_scores = []
        indirect_primers_scores = []
        # Intra strand: 
        for direct in direct_primers:
            primer = DNA_analyzer.reverseComp(self,strand[direct[0]-direct[1]:direct[0]+direct[1]]
                                      )
            intra53 = 1/np.mean([DNA_analyzer.measureComp(self,primer,strand[i-len(primer):i])
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            intra35 = 1/np.mean([DNA_analyzer.measureComp(self,primer,strand[i-len(primer):i][::-1])
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            inter53 = 1/np.mean([DNA_analyzer.measureComp(self,primer,
                      DNA_analyzer.reverseComp(self,strand[i-len(primer):i]))
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            inter35 = 1/np.mean([DNA_analyzer.measureComp(self,primer,
                      DNA_analyzer.reverseComp(self,strand[i-len(primer):i][::-1]))
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            direct_primers_scores.append(sum([intra53,intra35,inter53,inter35]))
            
        for indirect in indirect_primers:
            primer = strand[indirect[0]-indirect[1]:indirect[0]+indirect[1]]
            intra53 = 1/np.mean([DNA_analyzer.measureComp(self,primer,strand[i-len(primer):i])
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            intra35 = 1/np.mean([DNA_analyzer.measureComp(self,primer,strand[i-len(primer):i][::-1])
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            inter53 = 1/np.mean([DNA_analyzer.measureComp(self,primer,
                      DNA_analyzer.reverseComp(self,strand[i-len(primer):i]))
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            inter35 = 1/np.mean([DNA_analyzer.measureComp(self,primer,
                      DNA_analyzer.reverseComp(self,strand[i-len(primer):i][::-1]))
                      for i in range(len(primer),len(strand)-len(primer)-1)])
            indirect_primers_scores.append(sum([intra53,intra35,inter53,inter35]))
        direct = pd.DataFrame()
        direct["Start"] = [direct_primers[i][0]-direct_primers[i][1] for i in range(len(direct_primers))]
        direct["End"] = [direct_primers[i][0]+direct_primers[i][1] for i in range(len(direct_primers))]
        direct["Length"] = [direct_primers[i][1] for i in range(len(direct_primers))]
        direct["Specificity"] = [direct_primers_scores[i] for i in range(len(direct_primers_scores))]
        direct["Tm"] = [direct_primers[i][2] for i in range(len(direct_primers))]
        self.direct = direct
        indirect = pd.DataFrame()
        indirect["Start"] = [indirect_primers[i][0]-indirect_primers[i][1] for i in range(len(indirect_primers))]
        indirect["End"] = [indirect_primers[i][0]+indirect_primers[i][1] for i in range(len(indirect_primers))]
        indirect["Length"] = [indirect_primers[i][1] for i in range(len(indirect_primers))]
        indirect["Specificity"] = [indirect_primers_scores[i] for i in range(len(indirect_primers_scores))]
        indirect["Tm"] = [indirect_primers[i][2] for i in range(len(indirect_primers))]
        self.indirect = indirect
    
    # Show primers on the target sequence
    def showPrimers(self): 
        plt.style.use("dark_background")
        f = plt.figure(figsize=(15,10))
        ax = f.add_subplot(111)
        ax.plot([i for i in range(len(self.strand))],[1 for i in self.strand],c="blue")
        ax.plot([i for i in range(len(self.strand))],[-1 for i in self.strand],c="green")
        ax.set_ylim(-5,5)
        for i in range(self.direct.shape[0]):
            ax.hlines(xmin = self.direct.iloc[i,0],xmax = self.direct.iloc[i,1],y=3.5-(i+1)/10,color="red")
            if i==0:
                ax.annotate("".join(DNA_analyzer.reverseComp(self,self.strand[self.direct.iloc[i,0]:self.direct.iloc[i,1]])),(self.direct.iloc[0,1],4.5),color="red")
                ax.annotate("".join(self.strand[self.direct.iloc[i,0]:self.direct.iloc[i,1]]),(self.direct.iloc[0,1],4.2),color="blue")
        for i in range(self.indirect.shape[0]):
            ax.hlines(xmin = self.indirect.iloc[i,0],xmax = self.indirect.iloc[i,1],y=-3.5+(i+1)/10,color="red")
            if i==0:
                ax.annotate("".join(self.strand[self.indirect.iloc[i,0]:self.indirect.iloc[i,1]]),(self.indirect.iloc[0,1],-4.5),color="red")
                ax.annotate("".join(DNA_analyzer.reverseComp(self,self.strand[self.indirect.iloc[i,0]:self.indirect.iloc[i,1]])),(self.indirect.iloc[0,1],-4.2),color="green")
    def generate(self,length): 
        a = ["A","T","G","C"]
        return np.random.choice(a,length).tolist()
    
        
exp = DNA_primer_generator()
s = exp.generate(1000)
exp.FindPrimers(s,200,500)
exp.CheckPrimers()
exp.showPrimers()
