function result=inputnnn(s,t);
 result='';
 while isempty(result),
  if t=='num',
   result=input(s);
   else
   result=input(s,t);
   end;
  end;
