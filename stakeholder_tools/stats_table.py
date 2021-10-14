import numpy as np
import pandas as pd
import json
import panel as pn
import param
import os

from bokeh.models.widgets.tables import NumberFormatter, BooleanFormatter

class StatsTable():
    """StatsTable provides an interactive table to visualize the returned 
       statistics dictionary associated with `tclean` in the user stakeholder 
       tests.

       Args:
            json_measured (str): JSON dictionary associated with the measured 
            metric valuses from the `tclean` process.
            
            json_expected (str): JSON dictionary with the stakeholder defined 
            target metrics.
    """

    def __init__(self, json_measured:str, json_expected:str):
        pn.extension(sizing_mode="stretch_width")

        try:
            if os.path.isfile(json_measured) is False:
                bad_file = json_measured
                
                raise FileNotFoundError

            if os.path.isfile(json_expected) is False:
                bad_file = json_measured
                
                raise FileNotFoundError

        except FileNotFoundError:
            print('File doesn\'t exist: ' + bad_file)


            
        self.__json_measured = json_measured
        self.__json_expected = json_expected
        self.stakeholder_test = 'test_standard_cube'
        self.sub_test_stub = 'im_stats'
        self.df = None

    def __build_stats_dataframe(self)->pd.DataFrame:
        """Parses the stakeholder JSON files and builds a pandas DataFrame 
           with the reorganized metrics and keys. 

        Returns:
            pd.DataFrame: Pandas DataFrame containing organized keys, metrics
            from stakeholder tests.
        """

        keys = []
        expected = []
        measured = []
        
        with open(self.__json_measured) as file: 
            js = json.loads(file.read())
    

        for key, value in js[self.stakeholder_test][self.sub_test_stub + '_dict'].items():
            if isinstance(value, list):
                for v in value:
                    keys.append(key)
                    measured.append(v)
            else:
                keys.append(key)
                measured.append(value)
            
        with open(self.__json_expected) as file:
            js = json.loads(file.read())
    
        for key, value in js[self.stakeholder_test]['exp_' + self.sub_test_stub].items():
            if isinstance(value[-1], list):
                for v in value[-1]:
                    expected.append(v)
            else:
                expected.append(value[-1])
            
        self.df =  pd.DataFrame({
            'key': keys,
            'measured': measured,
            'expected': expected
        })

        return self.df

    @staticmethod
    def __color_negative_red(val)->str:
        """Takes a scalar and returns a string with
           the css property `'color: red'` for negative
           strings, black otherwise.

        Args:
            val ([float]): Pandas DataFrame values from the expected and measured columns. 

        Returns:
            [str]: string defining the value color properties.
        """
        
        color = 'gray' if val < 0 else 'black'
        
        return 'color: %s' % color
    
    @property
    def json_measured(self):
        """Getter funciton to get the measured metric values in JSON format.

        Returns:
            JSON object [JSON]: returns staekholder `measured` JSON object
        """
        with open(self.__json_measured) as file: 
            js = json.loads(file.read())

        return pn.pane.JSON(js, name='JSON', height=300, width=500)
    
    @json_measured.setter
    def json_measured(self, measured):
        """Setter funciton to set the measured metric file.
        """

        self.__json_measured = measured
    
    @property
    def json_expected(self):
        """Getter funciton to get the expected metric values in JSON format.

        Returns:
            JSON object [JSON]: returns staekholder `expected` JSON object
        """
        with open(self.__json_expected) as file: 
            js = json.loads(file.read())

        return pn.pane.JSON(js, name='JSON', height=300, width=500)

    @json_expected.setter
    def json_expected(self, expected):
        """Setter funciton to set the expected metric file.
        """

        self.__json_expected = expected

    def stats_table(self, tolerance = 0.05)->pn.widgets.Tabulator:
        """Builds and displays an interactive table containing the stakeholder tests results
           and expected metrics. The table contains two collapsible sections, one for passing
           values and one for failing values. 
        
        Args:
            tolerance ([float]): Tolerance for pass/fail of stakeholder test.

        Returns:
            pn.widgets.Tabulator: panel Tabulator widget
        """
    
        df = self.__build_stats_dataframe()
    
        df['pass'] = (abs((df.expected - df.measured)/df.expected)) < tolerance
    
        bokeh_formatters = {
            'value': NumberFormatter(format='0.00000'),
            'expected': NumberFormatter(format='0.00000'),
            'pass': BooleanFormatter(),
        }
    
        pn.widgets.Tabulator.theme = 'site'
    
        table = pn.widgets.Tabulator(
            df, 
            layout='fit_data_stretch', 
            formatters=bokeh_formatters,
            groupby=['pass'],
            pagination='remote'
        )
    
        return table

