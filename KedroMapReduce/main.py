from kedro.pipeline import Pipeline, node
from kedro.runner import SequentialRunner
from kedro.io import MemoryDataset, DataCatalog

from glob import glob
import json
import nltk
import re
import string


# nltk.download('punkt')


def get_words(raw_text: str) -> list[str]:
    words = []
    for tokens in nltk.sent_tokenize(raw_text):
        for word in nltk.word_tokenize(tokens):
            words.append(word.lower())
    return words


def map_words(raw_text: str) -> dict[str, int]:
    words = get_words(raw_text)
    result = {}
    for word in words:
        if word not in result:
            result[word] = 0
        result[word] += 1
    return result


def reduce_words(map_node1: dict[str, int], map_node2: dict[str, int], map_node3: dict[str, int], map_node4: dict[str, int], letters: list[str]) -> dict[str, tuple[str, int]]:
    result = {}
    all_keys = set(list(map_node1.keys()) + list(map_node2.keys()) + list(map_node3.keys()) + list(map_node4.keys()))
    for key in all_keys:
        if key[0] in letters:
            result[key] = map_node1.get(key, 0) + map_node2.get(key, 0) + map_node3.get(key, 0) + map_node4.get(key, 0)
    most_common_words = {}
    for letter in letters:
        result_letter = dict([(key, val) for key, val in result.items() if key[0].startswith(letter)])
        if not result_letter.keys():
            continue
        most_frequent_word = max(result_letter, key=result_letter.get)
        most_common_words[letter] = (most_frequent_word, result_letter[most_frequent_word])
    json.dump(most_common_words, open(''.join(letters) + '.json', 'w'), indent=1, ensure_ascii=False)
    return most_common_words

def main():
    files = sorted(glob('corpus/*'))
    all_text = ''
    for file in files:
        with open(file, 'r') as f:
            all_text += ' '
            all_text += ''.join(f.readlines()).strip()
    all_text = re.sub(f'[{string.punctuation}]', ' ', all_text)
    words = all_text.split()
    l = len(words)
    raw_text1 = ' '.join(words[:l//4])
    raw_text2 = ' '.join(words[l//4:l//2])
    raw_text3 = ' '.join(words[l//2:3*l//4])
    raw_text4 = ' '.join(words[3*l//4:])
        
    data_catalog = DataCatalog({
        'data1': MemoryDataset(raw_text1),
        'data2': MemoryDataset(raw_text2),
        'data3': MemoryDataset(raw_text3),
        'data4': MemoryDataset(raw_text4),
        'letters_ab': MemoryDataset(['a', 'b']),
        'letters_cd': MemoryDataset(['c', 'd']),
    })
    map_node1 = node(func=map_words, inputs=['data1'], outputs="map_node1")
    map_node2 = node(func=map_words, inputs=['data2'], outputs="map_node2")
    map_node3 = node(func=map_words, inputs=['data3'], outputs="map_node3")
    map_node4 = node(func=map_words, inputs=['data4'], outputs="map_node4")
    
    reduce_node1 = node(func=reduce_words, inputs=["map_node1", "map_node2", "map_node3", "map_node4", 'letters_ab'], outputs="reduce_node1")
    reduce_node2 = node(func=reduce_words, inputs=["map_node1", "map_node2", "map_node3", "map_node4", 'letters_cd'], outputs="reduce_node2")
    
    pipeline = Pipeline([map_node1, map_node2, map_node3, map_node4, reduce_node1, reduce_node2])

    runner = SequentialRunner()
    runner.run(pipeline, data_catalog)


if __name__ == '__main__':
    main()