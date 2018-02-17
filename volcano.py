from scipy import stats
import matplotlib   
import matplotlib.pyplot as plt

class volcano(object):

    def __init__(self, treated, control, expression_array=[], confidence_area=(None,0,'k',0), correction_type='FDR', color=('same', 'k'), size=('same',0,0)):
        self.treated = np.nan_to_num(treated)
        self.control = np.nan_to_num(control)
        self.expression_array = np.nan_to_num(expression_array)
        self.confidence_area = confidence_area
        self.correction_type = correction_type
        self.color = color
        self.size = size
        self.p_values = []          # to plot the p value, DO:  -1*np.log10(p)
        self.log2ratios = []
        self.fig = None
    
        # Check there are not Nans nor infs
        if expression_array==[]:
            self.expression_array = np.zeros(self.control.shape[0])
        all_arr = np.hstack([self.treated, self.control, self.expression_array.reshape(-1,1)])
        if True in np.isnan(all_arr) or True in np.isinf(all_arr):
            print('NaNs and infinite values found in your samples...')
        del all_arr
    

    def get_color(self):
        if self.color[0]=='same':
            cols = self.color[1]
        elif self.color[0]== 'expression':
            try:
                # In order to be fair in coloring, which will highly depend in the distribution of the data, I use here 
                # the 25% and 75% as min an max, so 50% of the data will be saturated in the lower or higher quintiles.
                (vmin, vmax) = ( np.array([0.20, 0.80]) * self.expression_array.shape[0] ).astype(int)
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
                cmap = eval('matplotlib.cm.' + self.color[1])
                cols = [cmap(norm(i)) for i in range(1,11)]

            except ValueError:
                print('color should be a tupple and don\' forget the expression array!')

        return cols



    def get_size(self):
        if self.size[0] == 'same':
            szes = 3
        else:
            szes = np.linspace(self.size[1].max(), self.size[1].min(), 10)
        
        return szes


    def get_p_values(self):
        # Check the shape of the arrays 
        shape_treated, shape_control = self.treated.shape, self.control.shape
        if shape_treated != shape_control:
            print('Treated and control arrays should have identical dimensions')
        spls = np.hstack([self.treated, self.control])

        # fil p_values and log2ratio arrays
        p_values, log2ratios = [], []
        for n_spl in range(len(self.control)):
            T, p = stats.ttest_ind(self.treated[n_spl], self.control[n_spl])
            ratio = np.mean(self.treated[n_spl]) / np.mean(self.control[n_spl])
            p_values.append(p)
            log2ratios.append(np.log2(ratio))

        self.p_values = np.array(p_values)
        self.log2ratios = np.array(log2ratios)


    def correction(self):
        if self.correction_type==None:
            pass
        elif self.correction_type=='FDR':
            pass
        elif self.correction_type=='FEWR':
            pass

    def show_confidence(self):
        if self.confidence_area[0]==None:
            return [0,0],[0,0]
        else:
            x_confidence, y_confidence = np.log2(self.confidence_area[0]), -np.log10(self.confidence_area[1])
            xMax, yMax = (-1 * log10(self.p_values)).max(), self.log2ratios.max()
            xMin, yMin = (-1 * log10(self.p_values)).min(), self.log2ratios.min()
        return [xMin*1.1, x_confidence, x_confidence, xMin*1.1], \
                [y_confidence, y_confidence, yMax*1.1, yMax*1.1] # this are xbox and ybox

    def PLOT(self):
        self.get_p_values()
        self.correction()
        
        colors = self.get_color()
        sizes = self.get_size()
        xbox, ybox = self.show_confidence()

        self.fig = plt.figure(figsize=(10,10))
        plt.scatter(self.log2ratios, -1*np.log10(self.p_values), s=sizes, c=colors) ## colorear por expression
        plt.fill(xbox, ybox, c=self.confidence_area[2], alpha=self.confidence_area[3])

        return self.fig

    @classmethod
    def plot(cls, treated, control, expression=[], relevant_areas=(None,0,'k',0), p_correction='FDR', color=('same', 'k'), size=('same',0,0)):
        volcano = cls(treated, control, expression, relevant_areas, p_correction, color, size)
        volcano.PLOT()
        plt.show()

        return volcano

    @staticmethod
    def help():
        print('Arguments: \n\
    Control array\n\
        - pandas or numpy, n-dimensions\n\
    Treated array\n\
        - pandas or numpy, n-dimensions\n\
    Expression_array\n\
    type of correction to the p-values: default is Benjamini-Hoschberg (\'FDR\')\n\
        - None \n\
        - Bonferroni (\'FWER\')\n\
        - Benjamini- Hochberg (\'FDR\')\n\
    color: Tuple. default is same color for all samples ((\'same\',\'hot\'))\n\
        - Control array should be already sorted by expression level (\'expression\')\n\
    size: Tuple. Default (\'same\', 0.1, 4)\n\
        - Control array should be already sorted by expression level (\'expression\')\n\
    confidence_area: default is None\n\
        - fill rectangle with level of confidence (tuple: (\'fold change\', p_value, color, alpha))\n')
