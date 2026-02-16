"""Code to parse file names for the ERA5 Land to DATM7 conversion scripts.

Functions
---------
resolve_file_sequences_or_format:
    Takes a sequence of file paths, or a callable, or a format string, and
    resolves it to a sequence of file paths, optionally checking for existence
    or non-existence, and/or duplicate values.
"""
from collections import Counter
from collections.abc import (
    Callable,
    Iterable,
    Mapping,
    Sequence,
)
from pathlib import Path
import typing as tp



def resolve_file_paths(
    paths: Sequence[Path] | str | Callable[..., Path],
    *,
    field_values: Mapping[str, Sequence[tp.Any]] | None = None,
    check_exists: bool = False,
    check_not_exists: bool = False,
    check_duplicates: bool = False,
) -> list[Path]:
    """
    Resolve file paths from various input formats.

    This function takes either a sequence of Path instances, a format string,
    or a callable, and returns a list of Path instances. For format strings
    and callables, keyword arguments are used to generate multiple paths by
    iterating over sequences of values.

    Parameters
    ----------
    paths : Sequence[Path] | str | Callable[..., Path]
        Either a sequence of Path instances, a format string with fields
        enclosed by `{}`, or a callable that returns a Path. NB! If a callable
        is provided, it should raise a `TypeError` if and only if it is called
        with missing or unexpected keyword arguments. Any TypeError raised will
        be interepreted as an indication that the keys in `field_values` do not
        match the parameters of the callable, and reported as such.
    field_values : Mapping[str, Sequence[tp.Any]] | None, optional
        Keyword arguments containing sequences of values. For format strings,
        these are used with `.format()`. For callables, these are passed as
        keyword arguments. All sequences must have the same length.
    check_exists : bool, optional
        If True, verify that all resolved paths exist. Default is False.
    check_not_exists : bool, optional
        If True, verify that none of the resolved paths exist. Default is False.
        No separate errors or warnings are raised if both check_exists and
        check_not_exists are True, and both checks will be performed, but at
        least one of them will of course fail.
    check_duplicates : bool, optional
        If True, verify that there are no duplicate paths in the resolved list.
        A DuplicateFilesError will be raised if any duplicates are found, with a
        list of the duplicate paths and their counts. Default is False.

    Returns
    -------
    list[Path]
        A list of resolved Path instances.

    Raises
    ------
    FilesNotFoundError
        If check_exists is True and any path does not exist, or if
        check_not_exists is True and any path exists, or if format_kwargs
        sequences have different lengths. NB! The raised exception is a
        FilesNotFoundError, which is a subclass of FileNotFoundError so that it
        will be caught by any code that catches FileNotFoundError, but it is a
        custom subclass with custom attributes.
    FilesAlreadyExistError
        If check_not_exists is True and any path exists. NB! The raised
        exception is a FilesAlreadyExistError, which is a subclass of
        FileExistsError so that it will be caught by any code that catches
        FileExistsError, but it is a custom subclass with custom attributes.
    DuplicateFilesError
        If check_duplicates is True and any duplicate paths are found.
    TypeError
        If paths is a format string or callable but `field_values` is not
        provided, or if the field names or keyword arguments expected by `paths`
        do not match the keys in `field_values`.

    Examples
    --------
    >>> # Using a sequence of paths
    >>> paths = resolve_file_paths([Path('file1.nc'), Path('file2.nc')])

    >>> # Using a format string
    >>> paths = resolve_file_paths(
    ...     'data_{year}_{month:02d}.nc',
    ...     year=[2020, 2020, 2021],
    ...     month=[1, 2, 1]
    ... )

    >>> # Using a callable
    >>> def make_path(year: int, month: int) -> Path:
    ...     return Path(f'data_{year}_{month:02d}.nc')
    >>> paths = resolve_file_paths(
    ...     make_path,
    ...     year=[2020, 2020, 2021],
    ...     month=[1, 2, 1]
    ... )
    """
    # Case 1: Sequence of Path instances
    if isinstance(paths, Sequence) and not isinstance(paths, str):
        resolved_paths: list[Path] = list(paths)
        missing_paths: set[Path]
        existing_paths: set[Path]

        if check_exists and not all(
                missing_paths := {
                    _path for _path in resolved_paths if not _path.exists()
                }
        ):
            raise FilesNotFoundError(missing_files=sorted(missing_paths))

        if check_not_exists and any(
                existing_paths := {
                    _path for _path in resolved_paths if _path.exists()
                }
        ):
            raise FilesAlreadyExistError(existing_files=sorted(existing_paths))

        return resolved_paths

    # Case 2 & 3: Format string or callable
    if field_values is None:
        raise TypeError(
            'When paths is a format string or callable, '
            'field_values must be provided'
        )
    if not isinstance(field_values, Mapping):
        raise TypeError(
            'field_values must be a mapping of field names to  field value '
            'sequences.'
        )

    # Verify all sequences have the same length
    lengths: set[int] = {len(_seq) for _seq in field_values.values()}
    if len(lengths) > 1:
        raise ValueError(
            f'All format_kwargs sequences must have the same length. '
            f'Got lengths::\n'
            + '\n'.join(
                f'- {_key}: {len(_seq)}'
                for _key, _seq in field_values.items()
            )
        )

    n_paths: int = lengths.pop() if lengths else 0

    # Generate paths by iterating over format_kwargs
    paths_func: Callable[..., Path]
    if isinstance(paths, str):
        paths_func: Callable[..., Path] = (
            lambda **kwargs: Path(paths.format(**kwargs))
        )
    else:
        paths_func: Callable[..., Path] = paths

    resolved_paths: list[Path] = [
        paths_func(**{_key: _seq[i] for _key, _seq in field_values.items()})
        for i in range(n_paths)
    ]

    # Perform optional checks
    if check_exists and not all(p.exists() for p in resolved_paths):
        missing: set[Path] = set(
            _path for _path in resolved_paths if not _path.exists()
        )
        raise FilesNotFoundError(missing_files=sorted(missing))

    if check_not_exists and any(p.exists() for p in resolved_paths):
        existing: set[Path] = set(
            _path for _path in resolved_paths if _path.exists()
        )
        raise FilesAlreadyExistError(existing_files=sorted(existing))

    if check_duplicates:
        path_counts: Counter[Path] = Counter(resolved_paths)
        duplicate_counts: dict[Path, int] = {
            _path: _count for _path, _count in path_counts.items() if _count > 1
        }
        if len(duplicate_counts) > 0:
            raise DuplicateFilesError(file_counts=duplicate_counts)

    return resolved_paths

###END def resolve_file_paths


class FilesNotFoundError(FileNotFoundError):
    """Custom error to indicate that no matching files were found

    Attributes
    ----------
    missing_files : list[Path]
        List of file paths that were expected but not found.
    """

    def __init__(self, *args, missing_files: Sequence[Path], **kwargs) -> None:
        """
        Parameters
        ----------
        *args
            Positional arguments to pass to the base FileNotFoundError.
        missing_files : Sequence[Path]
            List of file paths that were expected but not found.
        **kwargs
            Keyword arguments to pass to the base FileNotFoundError.
        """
        if len(args) == 0:
                args = (
                    (
                    f'No matching files were found. Missing files:\n'
                    + '\n'.join(f'- {file}' for file in missing_files)
                ),
            )
        super().__init__(*args, **kwargs)
        self.missing_files: tp.Final[list[Path]] = list(missing_files)
    ###END def __init__

###END class FilesNotFoundError


class FilesAlreadyExistError(FileExistsError):
    """Custom error to indicate that some files already exist

    Attributes
    ----------
    existing_files : list[Path]
        List of file paths that were expected not to exist but were found.
    """

    def __init__(self, *args, existing_files: Sequence[Path], **kwargs) -> None:
        """
        Parameters
        ----------
        *args
            Positional arguments to pass to the base FileExistsError.
        existing_files : Sequence[Path]
            List of file paths that were expected not to exist but were found.
        **kwargs
            Keyword arguments to pass to the base FileExistsError.
        """
        if len(args) == 0:
                args = (
                    (
                    f'The following files already exist:\n'
                    + '\n'.join(f'- {file}' for file in existing_files)
                ),
            )
        super().__init__(*args, **kwargs)
        self.existing_files: tp.Final[list[Path]] = list(existing_files)
    ###END def __init__

###END class FilesAlreadyExistError


class DuplicateFilesError(ValueError):
    """Custom error to indicate that duplicate file paths were found

    Attributes
    ----------
    duplicate_files : list[Path]
        List of file paths that were duplicated. If `file_counts` was passed to
        the constructor, this will be a list of all keys in that mapping that
        had a count greater than 1.
    file_counts : dict[Path, int] | None
        Optional dictionary mapping file paths to their counts in the input.
        This will be equal to a dictionary copied from the `file_counts`
        parameter in the constructor if it was provided, and None otherwise.
    """

    def __init__(
            self,
            *args,
            duplicate_files: Sequence[Path] | None = None,
            file_counts: Mapping[Path, int] | None = None,
            **kwargs
    ) -> None:
        """
        Parameters
        ----------
        *args
            Positional arguments to pass to the base ValueError.
        duplicate_files : Sequence[Path] | None, optional
            List of file paths that were duplicated. At least one of
            `duplicate_files` and `file_counts` must be provided. If
            `file_counts` is provided and is not None, then `duplicate_files`
            will be ignored, and the `duplicate_files` attribute will be set to
            a list of the keys in `file_counts` that had a count greater than 1.
        file_counts : Mapping[Path, int] | None, optional
            Dictionary mapping file paths to their counts in the input. At least
            one of `duplicate_files` and `file_counts` must be provided.
            `file_counts` overrides `duplicate_files` if both are provided.
        **kwargs
            Keyword arguments to pass to the base ValueError.
        """
        if duplicate_files is None:
            if file_counts is None:
                raise ValueError(
                    'At least one of duplicate_files and file_counts must be provided'
                )
            else:
                duplicate_files = [
                    _file for _file, _count in file_counts.items() if _count > 1
                ]
        if len(args) == 0:
            if file_counts is not None:
                counts_str = '\n'.join(
                    f'- {_file}: {_count}'
                    for _file, _count in file_counts.items() if _count > 1
                )
                args = (
                    (
                        f'The following files were duplicated, with counts:\n'
                        + counts_str
                    ),
                )
            else:
                args = (
                    (
                        f'The following files were duplicated:\n'
                        + '\n'.join(f'- {_file}' for _file in duplicate_files)
                    ),
                )
        super().__init__(*args, **kwargs)
        self.duplicate_files: tp.Final[list[Path]] = list(duplicate_files)
        self.file_counts: tp.Final[dict[Path, int] | None] = (
            dict(file_counts) if file_counts is not None else None
        )
    ###END def __init__

###END class DuplicateFilesError

