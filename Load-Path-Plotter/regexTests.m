endElements = "    NODE       SX             SY             SZ             SXY            SYZ            SXZ     ";

pat = "\s*(?<sectionKey>NODE)";

[match, nonMatch] = regexp(endElements, pat, 'names', 'split')
retval = true;
if ~isempty(match)
     if strcmp(match.sectionKey, "elem")
         if strcmp(match.context, "end")
             retVal = false;
         end
     end
end

