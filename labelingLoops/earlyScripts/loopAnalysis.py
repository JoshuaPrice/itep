import pandas as pd
import sqlite3 as sql3
import matplotlib.pyplot as plt
from epPlots import autolabel

def plotContactDistances():
	# Import df from sqlite
	sqlite_file = 'my_first_db.sqlite'    # name of the sqlite database file
	conn = sql3.connect(sqlite_file)
	c = conn.cursor()
	contactsDf = pd.read_sql_query("select * from AutoEP;", conn)
	conn.commit()
	c.close()
	conn.close()

	contactDistDf = contactsDf["startbpB"].subtract(contactsDf["endbpA"])

	fig, ax = plt.subplots()
	ax.set_xlabel('Contact Interaction Distance')
	ax.set_ylabel('Occurrence Count')
	ax.set_title('Distribution of contact distances in Rao 2014 CH12 cells')

	contactDistDf.plot.hist(bins=24, range=(0,1.2e6)) # I think spread is too wide for this - need to do more manually

	plt.show()


def plotTadJumps():
	# Import df from sqlite
	sqlite_file = 'my_first_db.sqlite'    # name of the sqlite database file
	conn = sql3.connect(sqlite_file)
	c = conn.cursor()
	contactsDf = pd.read_sql_query("select * from AutoEP;", conn)
	conn.commit()
	c.close()
	conn.close()

	tadsDf = contactsDf["tadA"]
	print tadsDf




if __name__ == "__main__":
	# recall TAD and reg labeled contacts file
	# plotContactDistances()
	plotTadJumps()


