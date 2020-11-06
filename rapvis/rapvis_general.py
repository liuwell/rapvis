#!/usr/bin/env python3
import datetime


###
def current_time():
	'''
	get current time
	'''
	return datetime.datetime.now().strftime('%b-%d-%Y %H:%M:%S')
