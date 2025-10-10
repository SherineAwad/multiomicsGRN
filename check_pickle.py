#!/usr/bin/env python3
import pickle
import sys
import os
from pprint import pprint

def debug_pickle(pkl_file):
    """
    Print everything possible from a pickle file for debugging
    """
    print(f"=== DEBUGGING PICKLE FILE: {pkl_file} ===")
    print(f"File exists: {os.path.exists(pkl_file)}")
    print(f"File size: {os.path.getsize(pkl_file) / 1024:.2f} KB")
    
    if not os.path.exists(pkl_file):
        print("ERROR: File does not exist!")
        return
    
    try:
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
    except Exception as e:
        print(f"ERROR loading pickle: {e}")
        return
    
    print(f"\n=== DATA TYPE: {type(data)} ===")
    
    if hasattr(data, '__len__'):
        print(f"Length: {len(data)}")
    
    # Print basic info based on data type
    if isinstance(data, dict):
        print_dict_details(data)
    elif isinstance(data, list):
        print_list_details(data)
    elif isinstance(data, tuple):
        print_tuple_details(data)
    elif hasattr(data, '__dict__'):  # Custom object
        print_object_details(data)
    else:
        print_generic_details(data)
    
    print(f"\n=== FULL CONTENTS ===")
    pprint(data, depth=3, width=120)

def print_dict_details(data):
    """Print details for dictionary"""
    print(f"Dictionary with {len(data)} keys:")
    for i, (key, value) in enumerate(data.items()):
        print(f"  [{i}] Key: {key} (type: {type(key)})")
        print(f"      Value: {type(value)} - {repr(value)[:100]}{'...' if len(repr(value)) > 100 else ''}")
        if i >= 10:  # Limit output for large dicts
            print(f"  ... and {len(data) - 10} more keys")
            break

def print_list_details(data):
    """Print details for list"""
    print(f"List with {len(data)} items:")
    for i, item in enumerate(data):
        print(f"  [{i}] {type(item)} - {repr(item)[:100]}{'...' if len(repr(item)) > 100 else ''}")
        if i >= 10:  # Limit output for large lists
            print(f"  ... and {len(data) - 10} more items")
            break

def print_tuple_details(data):
    """Print details for tuple"""
    print(f"Tuple with {len(data)} items:")
    for i, item in enumerate(data):
        print(f"  [{i}] {type(item)} - {repr(item)[:100]}{'...' if len(repr(item)) > 100 else ''}")

def print_object_details(data):
    """Print details for custom object"""
    print(f"Object of class: {data.__class__.__name__}")
    print(f"Module: {data.__class__.__module__}")
    print("Attributes:")
    for attr_name in dir(data):
        if not attr_name.startswith('__'):
            try:
                attr_value = getattr(data, attr_name)
                print(f"  {attr_name}: {type(attr_value)} - {repr(attr_value)[:80]}{'...' if len(repr(attr_value)) > 80 else ''}")
            except Exception as e:
                print(f"  {attr_name}: <cannot access - {e}>")

def print_generic_details(data):
    """Print details for generic data types"""
    print(f"Value: {repr(data)}")
    print(f"String representation: {str(data)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python debug_pickle.py <path_to_pickle_file>")
        sys.exit(1)
    
    pkl_file = sys.argv[1]
    debug_pickle(pkl_file)
