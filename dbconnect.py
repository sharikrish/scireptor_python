import configparser
import MySQLdb.cursors

config = configparser.ConfigParser()
config.read('/home/harikris-local/.my.cnf')

def connect():
    return MySQLdb.connect(host = config['mysql_igdb']['host'],
                           user = config['mysql_igdb']['user'],
      			   db = config['mysql_igdb']['database'])
