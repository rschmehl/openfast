"""
This was autogenerated with https://github.com/trolldbois/ctypeslib
and the header file generated from the registry.
I modified the original output to include some additional
tools for generating the C API in MoorDyn.
"""

import ctypes

# if local wordsize is same as target, keep ctypes pointer function.
if ctypes.sizeof(ctypes.c_void_p) == 8:
    POINTER_T = ctypes.POINTER
else:
    # required to access _ctypes
    import _ctypes
    # Emulate a pointer class using the approriate c_int32/c_int64 type
    # The new class should have :
    # ['__module__', 'from_param', '_type_', '__dict__', '__weakref__', '__doc__']
    # but the class should be submitted to a unique instance for each base type
    # to that if A == B, POINTER_T(A) == POINTER_T(B)
    ctypes._pointer_t_type_cache = {}
    def POINTER_T(pointee):
        # a pointer should have the same length as LONG
        fake_ptr_base_type = ctypes.c_uint64 
        # specific case for c_void_p
        if pointee is None: # VOID pointer type. c_void_p.
            pointee = type(None) # ctypes.c_void_p # ctypes.c_ulong
            clsname = 'c_void'
        else:
            clsname = pointee.__name__
        if clsname in ctypes._pointer_t_type_cache:
            return ctypes._pointer_t_type_cache[clsname]
        # make template
        class _T(_ctypes._SimpleCData,):
            _type_ = 'L'
            _subtype_ = pointee
            def _sub_addr_(self):
                return self.value
            def __repr__(self):
                return '%s(%d)'%(clsname, self.value)
            def contents(self):
                raise TypeError('This is not a ctypes pointer.')
            def __init__(self, **args):
                raise TypeError('This is not a ctypes pointer. It is not instanciable.')
        _class = type('LP_%d_%s'%(8, clsname), (_T,),{}) 
        ctypes._pointer_t_type_cache[clsname] = _class
        return _class

c_int128 = ctypes.c_ubyte*16
c_uint128 = c_int128
void = None
if ctypes.sizeof(ctypes.c_longdouble) == 16:
    c_long_double_t = ctypes.c_longdouble
else:
    c_long_double_t = ctypes.c_ubyte*16


"""
# class ProgDesc(Structure):
#     _fields_ = [
#         ("Name", POINTER(c_char)),
#         ("Ver", POINTER(c_char)),
#         ("Date", POINTER(c_char))
#     ]

Pointers must be listed first.

The fields here must be listed in the same order as the attributes in the
Fortran types. Basically, the way its listed in Fortran is the way it will
be used.

Whats with the padding in the autogenerated types?

Args:
    Structure ([type]): [description]
"""

class BaseOpenFASTType():
    def export_f2c(self):
        print( f"! {type(self).__name__} - Fortran to C" )
        for field in self._fields_:
            if "PADDING" in field[0]:
                continue
            print(f"{self._openfast_name}_C%{field[0]} = {self._openfast_name}%{field[0]}")

    def export_c2f(self):
        print( f"! {type(self).__name__} - C to Fortran" )
        for field in self._fields_:
            if "PADDING" in field[0]:
                continue
            print(f"{self._openfast_name}%{field[0]} = {self._openfast_name}_C%{field[0]}")
    
    def export_conversions(self):
        self.export_c2f()
        self.export_f2c()

class struct_MD_InitInputType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "InitInp"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('g', ctypes.c_float),
        ('rhoW', ctypes.c_float),
        ('WtrDepth', ctypes.c_float),
        ('PtfmInit', ctypes.c_float * 6),
        ('FileName', ctypes.c_char * 1024),
        ('RootName', ctypes.c_char * 1024),
        ('Echo', ctypes.c_bool),
        ('PADDING_0', ctypes.c_ubyte * 3),
        ('DTIC', ctypes.c_float),
        ('TMaxIC', ctypes.c_float),
        ('CdScaleIC', ctypes.c_float),
        ('threshIC', ctypes.c_float),
        ('OutList', POINTER_T(ctypes.c_char)),
        ('OutList_Len', ctypes.c_int32),
        ('PADDING_1', ctypes.c_ubyte * 4),
    ]

MD_InitInputType_t = struct_MD_InitInputType

# class struct_MD_LineProp(ctypes.Structure, BaseOpenFASTType):
#     _pack_ = True # source:False
#     _fields_ = [
#         ('object', POINTER_T(None)),
#         ('IdNum', ctypes.c_int32),
#         ('name', ctypes.c_char * 10),
#         ('PADDING_0', ctypes.c_ubyte * 2),
#         ('d', ctypes.c_float),
#         ('w', ctypes.c_float),
#         ('EA', ctypes.c_float),
#         ('BA', ctypes.c_float),
#         ('Can', ctypes.c_float),
#         ('Cat', ctypes.c_float),
#         ('Cdn', ctypes.c_float),
#         ('Cdt', ctypes.c_float),
#     ]
# MD_LineProp_t = struct_MD_LineProp

# class struct_MD_Connect(ctypes.Structure, BaseOpenFASTType):
#     _pack_ = True # source:False
#     _fields_ = [
#         ('object', POINTER_T(None)),
#         ('IdNum', ctypes.c_int32),
#         ('type', ctypes.c_char * 10),
#         ('PADDING_0', ctypes.c_ubyte * 2),
#         ('TypeNum', ctypes.c_int32),
#         ('PADDING_1', ctypes.c_ubyte * 4),
#         ('AttachedFairs', POINTER_T(ctypes.c_int32)),
#         ('AttachedFairs_Len', ctypes.c_int32),
#         ('PADDING_2', ctypes.c_ubyte * 4),
#         ('AttachedAnchs', POINTER_T(ctypes.c_int32)),
#         ('AttachedAnchs_Len', ctypes.c_int32),
#         ('conX', ctypes.c_float),
#         ('conY', ctypes.c_float),
#         ('conZ', ctypes.c_float),
#         ('conM', ctypes.c_float),
#         ('conV', ctypes.c_float),
#         ('conFX', ctypes.c_float),
#         ('conFY', ctypes.c_float),
#         ('conFZ', ctypes.c_float),
#         ('conCa', ctypes.c_float),
#         ('conCdA', ctypes.c_float),
#         ('Ftot', ctypes.c_float * 3),
#         ('Mtot', ctypes.c_float * 3 * 3),
#         ('S', ctypes.c_float * 3 * 3),
#         ('r', ctypes.c_float * 3),
#         ('rd', ctypes.c_float * 3),
#     ]
# MD_Connect_t = struct_MD_Connect

# class struct_MD_Line(ctypes.Structure, BaseOpenFASTType):
#     _pack_ = True # source:False
#     _fields_ = [
#         ('object', POINTER_T(None)),
#         ('IdNum', ctypes.c_int32),
#         ('type', ctypes.c_char * 10),
#         ('PADDING_0', ctypes.c_ubyte * 2),
#         ('OutFlagList', ctypes.c_int32 * 20),
#         ('FairConnect', ctypes.c_int32),
#         ('AnchConnect', ctypes.c_int32),
#         ('PropsIdNum', ctypes.c_int32),
#         ('N', ctypes.c_int32),
#         ('UnstrLen', ctypes.c_float),
#         ('BA', ctypes.c_float),
#         ('r', POINTER_T(ctypes.c_float)),
#         ('r_Len', ctypes.c_int32),
#         ('PADDING_1', ctypes.c_ubyte * 4),
#         ('rd', POINTER_T(ctypes.c_float)),
#         ('rd_Len', ctypes.c_int32),
#         ('PADDING_2', ctypes.c_ubyte * 4),
#         ('q', POINTER_T(ctypes.c_float)),
#         ('q_Len', ctypes.c_int32),
#         ('PADDING_3', ctypes.c_ubyte * 4),
#         ('l', POINTER_T(ctypes.c_float)),
#         ('l_Len', ctypes.c_int32),
#         ('PADDING_4', ctypes.c_ubyte * 4),
#         ('lstr', POINTER_T(ctypes.c_float)),
#         ('lstr_Len', ctypes.c_int32),
#         ('PADDING_5', ctypes.c_ubyte * 4),
#         ('lstrd', POINTER_T(ctypes.c_float)),
#         ('lstrd_Len', ctypes.c_int32),
#         ('PADDING_6', ctypes.c_ubyte * 4),
#         ('V', POINTER_T(ctypes.c_float)),
#         ('V_Len', ctypes.c_int32),
#         ('PADDING_7', ctypes.c_ubyte * 4),
#         ('T', POINTER_T(ctypes.c_float)),
#         ('T_Len', ctypes.c_int32),
#         ('PADDING_8', ctypes.c_ubyte * 4),
#         ('Td', POINTER_T(ctypes.c_float)),
#         ('Td_Len', ctypes.c_int32),
#         ('PADDING_9', ctypes.c_ubyte * 4),
#         ('W', POINTER_T(ctypes.c_float)),
#         ('W_Len', ctypes.c_int32),
#         ('PADDING_10', ctypes.c_ubyte * 4),
#         ('Dp', POINTER_T(ctypes.c_float)),
#         ('Dp_Len', ctypes.c_int32),
#         ('PADDING_11', ctypes.c_ubyte * 4),
#         ('Dq', POINTER_T(ctypes.c_float)),
#         ('Dq_Len', ctypes.c_int32),
#         ('PADDING_12', ctypes.c_ubyte * 4),
#         ('Ap', POINTER_T(ctypes.c_float)),
#         ('Ap_Len', ctypes.c_int32),
#         ('PADDING_13', ctypes.c_ubyte * 4),
#         ('Aq', POINTER_T(ctypes.c_float)),
#         ('Aq_Len', ctypes.c_int32),
#         ('PADDING_14', ctypes.c_ubyte * 4),
#         ('B', POINTER_T(ctypes.c_float)),
#         ('B_Len', ctypes.c_int32),
#         ('PADDING_15', ctypes.c_ubyte * 4),
#         ('F', POINTER_T(ctypes.c_float)),
#         ('F_Len', ctypes.c_int32),
#         ('PADDING_16', ctypes.c_ubyte * 4),
#         ('S', POINTER_T(ctypes.c_float)),
#         ('S_Len', ctypes.c_int32),
#         ('PADDING_17', ctypes.c_ubyte * 4),
#         ('M', POINTER_T(ctypes.c_float)),
#         ('M_Len', ctypes.c_int32),
#         ('LineUnOut', ctypes.c_int32),
#         ('LineWrOutput', POINTER_T(ctypes.c_float)),
#         ('LineWrOutput_Len', ctypes.c_int32),
#         ('PADDING_18', ctypes.c_ubyte * 4),
#     ]
# MD_Line_t = struct_MD_Line

# class struct_MD_OutParmType(ctypes.Structure, BaseOpenFASTType):
#     _pack_ = True # source:False
#     _fields_ = [
#         ('object', POINTER_T(None)),
#         ('Name', ctypes.c_char * 10),
#         ('Units', ctypes.c_char * 10),
#         ('QType', ctypes.c_int32),
#         ('OType', ctypes.c_int32),
#         ('NodeID', ctypes.c_int32),
#         ('ObjID', ctypes.c_int32),
#         ('PADDING_0', ctypes.c_ubyte * 4),
#     ]
# MD_OutParmType_t = struct_MD_OutParmType

class struct_MD_InitOutputType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "InitOut"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('writeOutputHdr', POINTER_T(ctypes.c_char)),
        ('writeOutputHdr_Len', ctypes.c_int32),
        ('PADDING_0', ctypes.c_ubyte * 4),
        ('writeOutputUnt', POINTER_T(ctypes.c_char)),
        ('writeOutputUnt_Len', ctypes.c_int32),
        ('PADDING_1', ctypes.c_ubyte * 4),
    ]
MD_InitOutputType_t = struct_MD_InitOutputType

class struct_MD_ContinuousStateType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "x"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('states', POINTER_T(ctypes.c_float)),
        ('states_Len', ctypes.c_int32),
        ('PADDING_0', ctypes.c_ubyte * 4),
    ]
MD_ContinuousStateType_t = struct_MD_ContinuousStateType

class struct_MD_DiscreteStateType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "xd"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('dummy', ctypes.c_float),
        ('PADDING_0', ctypes.c_ubyte * 4),
    ]
MD_DiscreteStateType_t = struct_MD_DiscreteStateType

class struct_MD_ConstraintStateType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "z"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('dummy', ctypes.c_float),
        ('PADDING_0', ctypes.c_ubyte * 4),
    ]
MD_ConstraintStateType_t = struct_MD_ConstraintStateType

class struct_MD_OtherStateType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "other"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('dummy', ctypes.c_float),
        ('PADDING_0', ctypes.c_ubyte * 4),
    ]
MD_OtherStateType_t = struct_MD_OtherStateType

class struct_MD_MiscVarType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "m"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('FairIdList', POINTER_T(ctypes.c_int32)),
        ('FairIdList_Len', ctypes.c_int32),
        ('PADDING_0', ctypes.c_ubyte * 4),
        ('ConnIdList', POINTER_T(ctypes.c_int32)),
        ('ConnIdList_Len', ctypes.c_int32),
        ('PADDING_1', ctypes.c_ubyte * 4),
        ('LineStateIndList', POINTER_T(ctypes.c_int32)),
        ('LineStateIndList_Len', ctypes.c_int32),
        ('PADDING_2', ctypes.c_ubyte * 4),
        ('MDWrOutput', POINTER_T(ctypes.c_float)),
        ('MDWrOutput_Len', ctypes.c_int32),
        ('PADDING_3', ctypes.c_ubyte * 4),
    ]
MD_MiscVarType_t = struct_MD_MiscVarType

class struct_MD_ParameterType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "p"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('NTypes', ctypes.c_int32),
        ('NConnects', ctypes.c_int32),
        ('NFairs', ctypes.c_int32),
        ('NConns', ctypes.c_int32),
        ('NAnchs', ctypes.c_int32),
        ('NLines', ctypes.c_int32),
        ('g', ctypes.c_float),
        ('rhoW', ctypes.c_float),
        ('WtrDpth', ctypes.c_float),
        ('kBot', ctypes.c_float),
        ('cBot', ctypes.c_float),
        ('dtM0', ctypes.c_float),
        ('dtCoupling', ctypes.c_float),
        ('NumOuts', ctypes.c_int32),
        ('RootName', ctypes.c_char * 1024),
        ('Delim', ctypes.c_char * 1),
        ('PADDING_0', ctypes.c_ubyte * 3),
        ('MDUnOut', ctypes.c_int32),
    ]
MD_ParameterType_t = struct_MD_ParameterType

class struct_MD_InputType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "u"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
    ]
MD_InputType_t = struct_MD_InputType

class struct_MD_OutputType(ctypes.Structure, BaseOpenFASTType):
    _openfast_name = "y"
    _pack_ = True # source:False
    _fields_ = [
        ('object', POINTER_T(None)),
        ('WriteOutput', POINTER_T(ctypes.c_float)),
        ('WriteOutput_Len', ctypes.c_int32),
        ('PADDING_0', ctypes.c_ubyte * 4),
    ]
MD_OutputType_t = struct_MD_OutputType

class struct_MD_UserData(ctypes.Structure, BaseOpenFASTType):
    _pack_ = True # source:False
    _fields_ = [
        ('MD_InitInput', MD_InitInputType_t),
        ('MD_InitOutput', MD_InitOutputType_t),
        ('MD_ContState', MD_ContinuousStateType_t),
        ('MD_DiscState', MD_DiscreteStateType_t),
        ('MD_ConstrState', MD_ConstraintStateType_t),
        ('MD_OtherState', MD_OtherStateType_t),
        ('MD_Misc', MD_MiscVarType_t),
        ('MD_Param', MD_ParameterType_t),
        ('MD_Input', MD_InputType_t),
        ('MD_Output', MD_OutputType_t),
    ]
MD_t = struct_MD_UserData

__all__ = \
    ['MD_Connect_t', 'MD_ConstraintStateType_t',
    'MD_ContinuousStateType_t', 'MD_DiscreteStateType_t',
    'MD_InitInputType_t', 'MD_InitOutputType_t', 'MD_InputType_t',
    'MD_LineProp_t', 'MD_Line_t', 'MD_MiscVarType_t',
    'MD_OtherStateType_t', 'MD_OutParmType_t', 'MD_OutputType_t',
    'MD_ParameterType_t', 'MD_t', 'struct_MD_Connect',
    'struct_MD_ConstraintStateType', 'struct_MD_ContinuousStateType',
    'struct_MD_DiscreteStateType', 'struct_MD_InitInputType',
    'struct_MD_InitOutputType', 'struct_MD_InputType',
    'struct_MD_Line', 'struct_MD_LineProp', 'struct_MD_MiscVarType',
    'struct_MD_OtherStateType', 'struct_MD_OutParmType',
    'struct_MD_OutputType', 'struct_MD_ParameterType',
    'struct_MD_UserData']
