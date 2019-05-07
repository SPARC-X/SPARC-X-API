import os

for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name.endswith('.pot'):
            f = open(os.path.join(root,name), 'r')
            txt = f.read()
            txt = txt.replace('D','E')
            f.close()
            f = open(os.path.join(root,name),'w')
            f.write(txt)
            f.close()
            print(name)
