import decimal

class PvalueHandler:

    def __init__(self):
        # initialize decimal math
        decimal.setcontext(decimal.BasicContext)
        decimal.getcontext().prec = 9
        decimal.getcontext().Emin = -10000000
        decimal.getcontext().Emax = 10000000
        # values to remember
        self.minimum_value_str = "0"
        self.maximum_value_str = "1"
        self.minimum_value = decimal.Decimal(self.minimum_value_str)
        self.maximum_value = decimal.Decimal(self.maximum_value_str)
        self.field_name = "p-values"

    # returns a decimal representation of the input string, or else raises an exception
    def parse_string (self,value_as_string):
        d = decimal.Decimal(value_as_string)
        if d.is_nan():
            raise Exception (f"Unable to convert '{value_as_string}' to a decimal value")
        if d <= self.minimum_value:
            raise Exception (f"Error: {self.field_name} must be greater than {self.minimum_value_str}")
        if d > self.maximum_value:
            raise Exception(f"Error: {self.field_name} can not be greater than {self.maximum_value_str}")
        return d

    # returns a float representation of the negative log of the input decimal, or raises an exception
    def neg_log_of_decimal (self,p_value):
        # prevent negative 0 output
        if p_value == self.maximum_value:
            return 0
        # prevent Inf output
        if p_value == self.minimum_value:
            return 999
        return float( -p_value.log10() )

