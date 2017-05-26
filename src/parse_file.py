from __future__ import division
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas


class FileParser(object):

	def __init__(self, input_file):
		self._input = input_file


	def _read_file(self):
		"""
		Read a file and convert it into dictionary based on chromosomes
		chromatin: 0, 1, 2
		"""
		output = {}
		with open(self._input, "r") as input_file:
			for line in input_file:
				line = line.split()
				chrom = line[0]
				if chrom not in output.keys():
					output[chrom] = []
				start_pos = line[1]
				end_pos = line[2]
				output[chrom].append([int(start_pos), int(end_pos)])
		return output
		# with open(feature_data+"chromatin"+".pickle", 'wb') as handle:
		# 	pickle.dump(output, handle, protocol=2)

	def _count_percent(self, chromatin_profile, ideogram):
		self._chrom_dict = {}
		self._open_percent = {}
		with open(ideogram, "r") as ideo:
			ideo.readline()
			for line in ideo:
				line = line.split()
				self._chrom_dict[line[0]] = int(line[2])
		# print chromatin_profile.keys()
		for key in chromatin_profile.keys():
			if key not in self._chrom_dict.keys(): continue
			count = 0
			for i in chromatin_profile[key]:
				count += i[1]-i[0]
			self._open_percent[key] = count
		print self._open_percent.values()
		return (sum(self._open_percent.values())/sum(self._chrom_dict.values()))*100

	def _merge_overlap(self, dict):
		"""
		merge all the overlapped lists 
		:param dict: output from read_file
		:return: dictionary
		"""
		for key in dict.keys():
			merged = []
			prev = [0,0]
			for item in dict[key]:
				if item[0] <= prev[1]:
					if item[1] <= prev[1]:
						continue
					prev[1] = item[1]
				else:
					merged.append(item)
					prev = item
			dict[key] = merged
		return dict

	def _get_intersect(self, l1, l2):
		"""
		get intersect of l1 and l2
		:param l1: 
		:param l2: 
		:return: 
		"""
		return list(set(l1) & set(l2))




# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
def chromosome_collections(df, y_positions, height,  **kwargs):
	"""
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
	del_width = False
	if 'width' not in df.columns:
		del_width = True
		df['width'] = df['end'] - df['start']
	for chrom, group in df.groupby('chrom'):
		print chrom
		yrange = (y_positions[chrom], height)
		xranges = group[['start', 'width']].values
		yield BrokenBarHCollection(
			xranges, yrange, facecolors=group['colors'], **kwargs)
	if del_width:
		del df['width']


def plot_profile(ideo_file, chromatin_file):

	# Height of each ideogram
	chrom_height = 1

	# Spacing between consecutive ideograms
	chrom_spacing = 1

	# Height of the gene track. Should be smaller than `chrom_spacing` in order to
	# fit correctly
	gene_height = 0.4

	# Padding between the top of a gene track and its corresponding ideogram
	gene_padding = 0.1

	# Width, height (in inches)
	figsize = (6, 8)

	# Decide which chromosomes to use
	chromosome_list = ['chr%s' % i for i in range(1, 23) + ['M', 'X', 'Y']]

	# Keep track of the y positions for ideograms and genes for each chromosome,
	# and the center of each ideogram (which is where we'll put the ytick labels)
	ybase = 0
	chrom_ybase = {}
	gene_ybase = {}
	chrom_centers = {}

	# Iterate in reverse so that items in the beginning of `chromosome_list` will
	# appear at the top of the plot
	for chrom in chromosome_list[::-1]:
		chrom_ybase[chrom] = ybase
		chrom_centers[chrom] = ybase + chrom_height / 2.
		gene_ybase[chrom] = ybase - gene_height - gene_padding
		ybase += chrom_height + chrom_spacing

	# Read in ideogram.txt, downloaded from UCSC Table Browser
	ideo = pandas.read_table(ideo_file,
		skiprows=1,
		names=['chrom', 'start', 'end', 'name', 'gieStain']
	)

	# Filter out chromosomes not in our list
	ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

	# Add a new column for width
	ideo['width'] = ideo.end - ideo.start

	# Colors for different chromosome stains
	color_lookup = {
		'gneg': "lightgrey",
		'gpos25': (1., 1., 1.),
		'gpos50': "lightgrey",
		'gpos75': (1., 1., 1.),
		'gpos100': (1., 1., 1.),
		'acen': (1., 1., 1.),
		'gvar': "lightgrey",
		'stalk': (1., 1., 1.),
	}

	# Add a new column for colors
	ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])
	# ideo['colors'] = 'white'

	# Same thing for genes
	genes = pandas.read_table(
		chromatin_file,
		names=['chrom', 'start', 'end', "name"],
		usecols=range(4))
	genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
	genes['width'] = genes.end - genes.start
	genes['colors'] = "lightgrey"


	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)

	# Now all we have to do is call our function for the ideogram data...
	print("adding ideograms...")
	for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
		ax.add_collection(collection)

	# ...and the gene data
	print("adding genes...")

	for collection in chromosome_collections(genes, gene_ybase, gene_height, alpha=5, linewidths=0.1):
		ax.add_collection(collection)

	# Axes tweaking
	ax.set_yticks([chrom_centers[i] for i in chromosome_list])
	ax.set_yticklabels(chromosome_list)
	ax.axis('tight')
	plt.title("Breast")
	plt.show()


if __name__ == "__main__":
	ideo_file = '/Users/roujia/Documents/02_dev/09_chromatin_accessibility/src/cytoBandIdeo.txt'
	chromatin_file = '/Users/roujia/Documents/02_dev/09_chromatin_accessibility/Data/breast-ENCFF001UVV.bed'

	# plot one profile
	# plot_profile(ideo_file, chromatin_file)

	# test merge
	# test = {"a": [[1, 2], [8, 9], [10,16], [15,16]], "b": [[1, 2], [4,10],[19,20],[5,9]]}
	obj = FileParser("/Users/roujia/Documents/02_dev/09_chromatin_accessibility/Data/skin-fetal-ENCFF976JED.bed")
	# obj._merge_overlap(test)
	chromatin_profile = obj._read_file()
	chromatin_profile = obj._merge_overlap(chromatin_profile)
	print obj._count_percent(chromatin_profile, ideo_file)
