import pandas as pd

path = r"C:\Users\Ivan\Desktop\THESIS\Sources\workhours.csv"

df = pd.read_csv(path)

print(df.head())