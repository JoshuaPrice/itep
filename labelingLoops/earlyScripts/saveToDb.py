import pandas as pd 
import sqlite3 as sql3

sqlite_file = 'my_first_db.sqlite'    # name of the sqlite database file

# Connecting to the database file
conn = sql3.connect(sqlite_file)
c = conn.cursor()

# Load df
labeledLoopsFile = '/data2/josh/expCH12/automated_labeled_loops.csv'
labeledLoopsDf = pd.read_csv(labeledLoopsFile, sep="\t")

# Save df into sqlite table
labeledLoopsDf.to_sql("AutoEP", conn, if_exists="replace")
conn.commit()

#Query database
print pd.read_sql_query("select * from AutoEP;", conn)

# Close cursor and connection
c.close()
conn.close()


