import json

file = open('overall.json', 'r')
overall_list = json.loads(file.read())
print(len(overall_list))
print(overall_list[0]['molecules'][0]['location'])