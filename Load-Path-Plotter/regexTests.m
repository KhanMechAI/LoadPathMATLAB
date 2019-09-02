tesString = "   NODE       SX             SY             SZ             SXY            SYZ            SXZ     ";

pat = "\s*(?<sectionKey>NODE)\s+(?:(?<stressDim>S\w+)\s*)+";
pat = "[\s\t\n\r]*(NODE)[\s\t\n\r]*((S\w+)[\s\t\n\r]*)+";
pat = "\s*(?<sectionKey>NODE)\s*(?<stressFields>\<S(?<sub>\w+))*"
pat = ["\s+(?<sectionKey>NODE)\s+","(?<stressFields>\<S\w+)\s+"]
[match, nonMatch] = regexp(tesString, pat, 'names')
retval = true
if ~isempty(match)
    if length(match)>1
        stressFields = match(1,2);
        stressFields
    end
end

