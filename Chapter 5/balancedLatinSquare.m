% Based on "Bradley, J. V. Complete counterbalancing of immediate sequential effects in a Latin square design. J. Amer. Statist. Ass.,.1958, 53, 525-528. "
% How to use:
% var conditions = ["A", "B", "C", "D"];
% balancedLatinSquare(conditions, 0)  => ["A", "B", "D", "C"]
% balancedLatinSquare(conditions, 1)  => ["B", "C", "A", "D"]
% balancedLatinSquare(conditions, 2)  => ["C", "D", "B", "A"]
% ...
function result = balancedLatinSquare(array, participantID)

    result = [];
    num_of_conditions = length(array);
    j = 0;
    h = 0;
    
    for i = 0: num_of_conditions - 1
        val = 0;
        
        if(i < 2 || mod(i, 2) ~= 0)            
            val = j;
            j = j + 1;
        else
            val = num_of_conditions - h - 1;
            h = h + 1;
        end
        
        idx = mod((val + participantID), num_of_conditions) + 1;
        result = [result array(idx)];       
    end
    
    if(mod(num_of_conditions, 2) ~= 0 && mod(partipantID, 2) ~= 0)
        result = fliplr(result);
    end

end