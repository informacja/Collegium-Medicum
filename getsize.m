 exportgraphics(gcf, "a.pdf",...
'BackgroundColor','none', "contentType","vector")

sum = 0;
GetSize(gcf().CurrentObject)
for(i = 1:length(gcf().Children))
    
a = GetSize(gcf().Children(i));
sum = sum + a;
end
% Children
sum

function [totSize] = GetSize(this) 
   props = properties(this); 
   totSize = 0; 
   
   for ii=1:length(props) 
      currentProperty = getfield(this, char(props(ii))); 
      s = whos('currentProperty'); 
      totSize = totSize + s.bytes; 
   end
  
   % fprintf(1, '%d bytes\n', totSize);
   
end