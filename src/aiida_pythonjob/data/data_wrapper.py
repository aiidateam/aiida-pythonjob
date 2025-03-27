from aiida import orm


class Dict(orm.Dict):
    @property
    def value(self):
        return self.get_dict()


class List(orm.List):
    @property
    def value(self):
        return self.get_list()


class ArrayData(orm.ArrayData):
    @property
    def value(self):
        return self.get_array()
