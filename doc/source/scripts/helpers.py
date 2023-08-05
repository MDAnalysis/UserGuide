from typing import Optional, Callable, Type

def sphinx_class(*, klass: Type, tilde: bool =True) -> str:
    prefix = '~' if tilde else ''
    return ':class:`{}{}.{}`'.format(prefix, klass.__module__, klass.__name__)

def sphinx_meth(*, meth: Callable, tilde: bool=True) -> str:
    prefix = '~' if tilde else ''
    return ':meth:`{}{}.{}`'.format(prefix, meth.__module__, meth.__qualname__)

def sphinx_ref(*, txt: str, label: Optional[str] =None, suffix: str ='') -> str:
    if label is None:
        label = txt
    return ':ref:`{} <{}{}>`'.format(txt, label, suffix)

def sphinx_link(*, txt: str) -> str:
    return '`{}`_'.format(txt)
