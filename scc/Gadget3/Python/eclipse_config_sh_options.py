# usage: python eclipse_config_sh_options.py
#
# in the Eclipse menu: File -> Import -> C/C++ -> C/C++ Project Settings
# there choose eclipse_symbols.xml and import it
#
# to activate folding of undefined symbols in Eclipse
# in the Eclipse menu: Window -> Preferences -> C/C++ -> Editor -> Folding
# there click:
# "Enable folding when opening a new editor"
# "Enable folding of preprocessor branches"
# "Inactive preprocessor branches"
#
# to check if importing worked corectly: Project -> Properties -> C/C++ Build -> GCC C Compiler -> Symbols
 
import sys
import string

#text=sys.stdin.readlines() # usage in this case: cat ../Config.sh | python eclipse_config_sh_options.py > eclipse_symbols.xml

file_in = open("../Config.sh",'r')
text = file_in.readlines()

file_out = open("eclipse_symbols.xml",'w')

header = [
'<?xml version="1.0" encoding="UTF-8"?>',
'<cdtprojectproperties>',
'<section name="org.eclipse.cdt.internal.ui.wizards.settingswizards.Macros">',
'<language name="C Source File">',
'']
#sys.stdout.write (string.join(header, '\n'))
file_out.write(string.join(header, '\n'))

for curtok in text:
  curline = string.strip(curtok)
  if curline.__len__() >= 1:
    if curline[0]!="#":
      curline = (curline.split())[0]
      lines = ['<macro>',
               '<name>' + string.strip(curline) + '</name><value/>',
               '</macro>', '']
      #sys.stdout.write(string.join(lines, '\n'))
      file_out.write(string.join(lines, '\n'))

footer = [
'</language>',
'<language name="Object File">',
'',
'</language>',
'<language name="Assembly Source File">',
'',
'</language>',
'</section>',
'</cdtprojectproperties>',
'']
#sys.stdout.write (string.join(footer, '\n'))
file_out.write (string.join(footer, '\n'))
