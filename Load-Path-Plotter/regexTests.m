endElements = "    NODE       SX             SY             SZ             SXY            SYZ            SXZ     ";

pat = "\s+(?<sectionKey>NODE)\s+(?:(?<stressDim>S\w+)\s*)+";

[match, nonMatch] = regexp(endElements, pat, 'names', 'split')
retval = true;
if ~isempty(match)
     if strcmp(match.sectionKey, "elem")
         if strcmp(match.context, "end")
             retVal = false;
         end
     end
end

