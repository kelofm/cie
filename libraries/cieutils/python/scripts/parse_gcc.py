# --- STL Imports ---
import sys


class NamedObject:
    _invalidCharacters = " .+-*/%='\"!<>()"

    def __init__(self, name: str, validate=True):
        self._name, self._nameSpace = self._parseName(name, validate=validate)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name: str):
        self._name, self._nameSpace = self._parseName(name)

    def __str__(self):
        return "".join("{}::".format(nameSpace) for nameSpace in self._nameSpace) + self._name

    @property
    def fullName(self):
        return self.name

    @property
    def shortName(self):
        return self._name

    @staticmethod
    def _parseName(name: str, validate=True):
        """Split the input name to its namespace and name stub."""
        if validate:
            for char in NamedObject._invalidCharacters:
                if char in name:
                    raise NameError("Invalid character '{}' in name: {}".format(char, name))

        if not name:
            return "", []
        
        name = name.replace("typename","").split("::")
        return name[-1], name[:-1]


class Template(NamedObject):
    def __init__(self, string: str, validate=True):
        NamedObject.__init__(self, "", validate=False)
        self._parameters = [] # set of parameters as tuples (alias,Template) 

        # Validate
        string = string.strip()

        # Find scope
        parametersBegin = string.find('<')
        parametersEnd   = string.rfind('>')

        if string and parametersBegin!=-1 and parametersEnd!=-1 and parametersBegin!=parametersEnd:
            # Set name
            self.name = string[:parametersBegin]

            # Parse parameters
            parameters = string[parametersBegin+1:parametersEnd].split(',')
            for parameter in parameters:
                parameter = parameter.strip()
                self._parameters.append(("", Template(parameter)))

        else:
            # Set name
            self.name = string

    
    def __str__(self):
        if self._parameters:
            return "{name}<{parameters}>".format(
                name = NamedObject.__str__(self),
                parameters = ",".join(str(pair[1]) for pair in self._parameters)
            )
        else:
            return NamedObject.__str__(self)



def extractEnclosedString(string: str, beginCharacter='(', endCharacter=')', validate=True):
    """beginCharacter contents endCharacter"""
    begin  = string.find(beginCharacter)
    end    = string.rfind(endCharacter)

    if begin == -1 or end == -1 or end < begin:
        if validate:
            raise ValueError("'{}' is not enclosed by '{}' and '{}'".format(string, beginCharacter, endCharacter))
    else:
        return string[begin+1:end].strip()

    return ""



class TemplateDetails:
    def __init__(self, string: str, validate=True):
        """[with shorthand_0 = expansion_0; shorthand_1 = expansion_1; ...]"""
        self._types = {}
        string = extractEnclosedString(string, beginCharacter='[', endCharacter=']', validate=validate)

        # Remove 'with'
        begin = string.find("with")
        if begin == -1:
            if validate:
                raise ValueError("Invalid template expansion '{}'".format(string))
            begin = 0
        else:
            begin += len("with")
        
        string = string[begin:]

        # Parse aliases
        while string:
            # Find scope
            end = string.find(';')
            if end == -1:
                end = len(string)

            contents = string[:end].split('=')
            string   = string[end+1:].strip()

            # Parse components
            if len(contents) != 2:
                raise RuntimeError("Invalid template expansion {}".format(contents))

            self._types[Template(contents[0])] = Template(contents[1])


    def __getitem__(self, key: str):
        return self._types[key]


    def keys(self):
        return self._types.keys()


    def values(self):
        return self._types.values()


    def items(self):
        return self._types.items()


    def __iter__(self):
        return self._types.__iter__()


    def __next__(self):
        return self._types.__next__()



class ArgumentList:
    def __init__(self, string: str, validate=True):
        """(argument_0_type, argument_1_type, ...)"""
        string = extractEnclosedString(string, validate=validate)
        self._types = string.split(',')


    def __getitem__(self, index: int):
        return self._types[index]


    def __iter__(self):
        return self._types.__iter__()


    def __next__(self):
        return self._types.__next__()




class Function(Template):
    def __init__(self, string: str, validate=True):
        """function_name(...arguments) [template_details] (not_yet_recognized)"""
        Template.__init__(self, "", validate=False)
        self._arguments      = ArgumentList("", validate=False)
        self_templateDetails = TemplateDetails("", validate=False)
        self._TODO           = ""

        # Find component separators
        argumentsBegin = string.find('(')
        argumentsEnd   = string.find(')') + 1
        detailsBegin   = string.find('[')
        detailsEnd     = string.find(']') + 1
        returnBegin    = string.rfind('(') + 1
        returnEnd      = string.rfind(')')

        # Parse name
        if argumentsBegin == -1:
            raise ValueError("Missing function name in '{}'".format(string))
        self.name = string[:argumentsBegin]

        # Parse arguments
        if argumentsBegin == -1 or argumentsEnd == 0 or argumentsEnd < argumentsBegin:
            raise ValueError("Missing function arguments in '{}'".format(string))
        self._arguments = ArgumentList(string[argumentsBegin:argumentsEnd])

        # Parse template expansion [optional]
        self._templateDetails = TemplateDetails(string[detailsBegin:detailsEnd], validate=False)

        # Parse return type
        if returnBegin == -1 or returnEnd == 0 or returnEnd < returnBegin:
            raise ValueError("Missing return type in '{}'".format(string))
        self._TODO = string[returnBegin, returnEnd]




if __name__ == "__main__":
    contents = ""

    #filePath = sys.argv[1]
    filePath = "/home/mate/Desktop/error.log"

    with open(filePath, "r") as file:
        contents = file.read().strip()

    print(Function(contents))