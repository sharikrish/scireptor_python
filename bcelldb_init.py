#!/usr/bin/python3

"""
@Srilakshmy - s.harikrishnan@dkfz-heidelberg.de
@franasa - f.arcila@dkfz-heidelberg.de
    -adapted from Peter Arndt and Christian Busse-

DESCRIPTION

- logging processes to the database; every calling of the module logs the calling main method (user, machine, command line arguments, ...)
- returning database handles (with configurations in the .my.cnf file in your home directory)

"""

import re
from configobj import ConfigObj # install via pip3
import MySQLdb.cursors # install via pip3
from os.path import expanduser
import logging

config=(ConfigObj(expanduser('~/.my.cnf')))

# logging.basicConfig(filename='example.log',level=logging.DEBUG)#

locals()
def connect():
    return MySQLdb.connect(host = config['mysql_igdb']['host'],
                           user = config['mysql_igdb']['user'],
      			   db = config['mysql_igdb']['database'])


re_inline_comment = re.compile("^(.*?)(?<!\\\\)#.*") # Match everything before a non-escaped hash
re_key_value = re.compile("^\s*([_A-Za-z][_0-9A-Za-z]+)=(.*?)\s*;?\s*$")


# try to open config file in .

def get_config():
    """
    Look for config file in . and than ../
    Return ConfigObj dictionary.
    """
    # try to open config file in .
    try:
        config_file = open("config", "r")

    except IOError:
        # try from ../ directory
        try:
            config_file = open("../config", "r")
        except IOError:
            print("no config file found")

    return ConfigObj(config_file)


if __name__ == "__main__":
    connect()
    get_config()
