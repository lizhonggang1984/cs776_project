from DBSCAN.dbscan import dbscan

def runDBSCAN():
    dbscan('expr_matrix.txt')

if __name__=='__main__':
    runDBSCAN()