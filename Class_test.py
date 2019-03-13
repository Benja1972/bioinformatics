class One:
    def __init__(self,ai):
        self.a = ai
    
    @property
    def temp(self):
        return self.a
    
    @temp.setter
    def temp(self,b):
        self.a=b


on = One(3)


on.temp=10


print(on.temp)
