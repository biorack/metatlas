
******************
Information
******************

Data Model
----------
The raw data is stored in mzML files at NERSC.
The raw data is parsed into HDF files for easier access and stored
alongside the mzML files.
The metatlas metadata is stored in a MySQL database at NERSC.
There are python objects in metabolite atlas that correspond to each
object type in the table.
The python objects are IPython (traitlets)[http://traitlets.readthedocs.org/en/stable/],
which allow for attribute validation and change notification.
The (dataset)[https://dataset.readthedocs.org/en/latest/] library is used
where possible to simplify access to the MySQL database.
Otherwise, custom queries are created as needed.
The dataset library provides custom table schema creation and ORM handling.

An object holds a reference to another object by storing the `unique_id`
of the other object in the table.  An object can contain a list
of references to other objects.  These are stored in a separate link table, which
allows a many-to-many relationship of objects. Where possible, access to the database
is batched, where a query to each table is made only once for a given operation.

Provenance Model
----------------
A metatlas object has a `creation_time` and a `username`,
which is the `$USER` of the person who created this version of the object.
An object has three unique identifiers: the `unique_id`, the `prev_uid`, and the `head_id`.
The "current object" is the one where the `unique_id` matches the `head_id`.
When a user makes a change to their own object, the previous entry in the
database is given a new `unique_id`, that id is stored as the `prev_uid` of the new entry.
The same is done for any entries in the link table(s).  The `head_id` is used for grouping
a set of changes together.  A new `head_id` is created on a `.clone()` of an object,
or when saving changes to an object created by another user.

Visualization
-------------
A metatlas object can be edited by calling `.edit()`, which creates an IPython widget
interface to the object, allowing the user to edit the object.
A group of metatlas objects can be edited by calling `edit_objects()`, which creates
a QGrid widget, allowing the user to edit the objects in a table.
The matplotlib notebook backend is used to create interactive visualizations of
data.
