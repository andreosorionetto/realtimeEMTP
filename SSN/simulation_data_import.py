
def get_column(name,arq):
    text=open(arq ,'r')
    column=[]
    lines=text.read().split('\n')
    space_bars=[]
    for n in range(0,len(lines[1])):
        if lines[1][n] == ' ' and n == 0:
            space_bars.append(n)
            
        if lines[1][n] == ' ' and n!= 0 and lines[1][n-1] != ' ':
            space_bars.append(n)
    
    names = lines[0].split(' ') 
    while names.count('')!=0:
        names.remove('')
    
    Dic={}
    for n in names:
        Dic[n] = (space_bars[names.index(n)], space_bars[names.index(n)+1]) # O intervalo que é preciso pegar para se obter aquela coluna
    
    i = Dic[name][0]
    f = Dic[name][1]
    
    for line in lines:
        column.append(line[ i : f ])
    text.close()
    return column
   
def create_signal(L): #Transforma em lista de numeros
    L_new=[]
    while L.count(' ') != 0:
        L.remove(' ')
    for e in L:
        try:
            L_new.append(float(e))
        except:
            pass
    return L_new



