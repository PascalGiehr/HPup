'''
Created on 10.08.2016

@author: Andy
'''
import configparser
import sys



    
config= configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
# first argument is the config file, second argument is the section, third argument is the parameter
config.read(sys.argv[1])

#key1 is section, key2 is Parameter
#if key1 exist in config file, print it else throw ERROR
#if key2 exist in config[key1], print it else throw Error 
def getArgument():
    config= configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(sys.argv[1])
    if config.has_section(sys.argv[2]):
        if config.has_option(sys.argv[2], sys.argv[3]):
            print(config.get(sys.argv[2], sys.argv[3]))
        else:
            print("Searched option does not exist!")
            sys.exit(1)
    else:
        print("Searched section does not exist!")
        print(sys.argv[1])
        print(sys.argv[2])
        print(sys.argv[3])
        sys.exit(1)

if __name__ == "__main__":
    getArgument()
    
