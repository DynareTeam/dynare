function bool=isvar(name)
%Returns true if there is a variable named NAME in the current workspace.

bool=evalin('caller', ['exist(''',name,''',''var'');']);