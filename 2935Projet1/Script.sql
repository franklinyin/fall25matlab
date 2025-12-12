-- create schema if not existing
begin;
CREATE SCHEMA IF NOT EXISTS projet;
commit;


begin;


-- drop order is important
DROP TABLE IF EXISTS projet.Borrow;
DROP TABLE IF EXISTS projet.Order;
DROP TABLE IF EXISTS projet.Subscriber;
DROP TABLE IF EXISTS projet.Books;
DROP TABLE IF EXISTS projet.Address;
DROP TABLE IF EXISTS projet.Borrow_Date;
DROP TABLE IF EXISTS projet.Return_Date;
DROP TABLE IF EXISTS projet.Library;


-- table creation section
CREATE TABLE projet.Books (
    BookID SERIAL PRIMARY KEY,
    title VARCHAR(255),
    author VARCHAR(255),
    genre VARCHAR(255)
);


CREATE TABLE projet.Subscriber (
    SubscriberID SERIAL PRIMARY KEY,
    name VARCHAR(255),
    AddressID int,
    LibraryID int
);

CREATE TABLE projet.Address (
    AddressID SERIAL PRIMARY KEY,
    street VARCHAR(255),
    apartment int,
    postal_code VARCHAR(255),
    province VARCHAR(255),
    city VARCHAR(255)
);

CREATE TABLE projet.Borrow (
    BorrowID SERIAL PRIMARY KEY,
    Borrow_DateID int,
    Return_DateID int,
    duration int,
    possible_delay VARCHAR(255),
    BookID int,
    SubscriberID int,
    LibraryID int
);

CREATE TABLE projet.Borrow_Date (
    Borrow_DateID SERIAL PRIMARY KEY,
    day int,
    month int,
    year int
);

CREATE TABLE projet.Return_Date (
    Return_DateID SERIAL PRIMARY KEY,
    day int,
    month int,
    year int
);


CREATE TABLE projet.Order (
    BookID int,
    SubscriberID int,
    status VARCHAR(255),
    LibraryID int,
    PRIMARY KEY (BookID, SubscriberID)
);



CREATE TABLE projet.Library (
    LibraryID SERIAL PRIMARY KEY,
    name VARCHAR(255)
);



commit;


------------------------------------------------------------------------------------
begin;

-- data insertion section
INSERT INTO projet.Books (title, author, genre) VALUES
('Le Petit Prince'                                      , 'Antoine de Saint-Exupéry', 'Roman'),
('Harry Potter à l''école des sorciers'                 , 'J.K. Rowling'            , 'Fantasy'),
('1984'                                                 , 'George Orwell'           , 'Science-fiction'),
('Orgueil et Préjugés'                                  , 'Jane Austen'             , 'Roman classique'),
('Le Seigneur des Anneaux : La Communauté de l''Anneau' , 'J.R.R. Tolkien'          , 'Fantasy'),
('Crime et Châtiment'                                   , 'Fiodor Dostoïevski'      , 'Roman'),
('Les Misérables'                                       , 'Victor Hugo'             , 'Roman classique'),
('Le Vieil Homme et la Mer'                             , 'Ernest Hemingway'        , 'Roman'),
('L''Étranger'                                          , 'Albert Camus'            , 'Roman philosophique'),
('Le Hobbit'                                            , 'J.R.R. Tolkien'          , 'Fantasy');


INSERT INTO projet.Subscriber (name, AddressID, LibraryID) VALUES
('Jean Dupont'          , 1 , 1),
('Marie Tremblay'       , 2 , 1),
('Pierre Gagnon'        , 3 , 1),
('Sophie Martin'        , 4 , 1),
('François Leblanc'     , 5 , 1),
('Isabelle Bergeron'    , 6 , 1),
('Éric Lavoie'          , 7 , 1),
('Valérie Roy'          , 8 , 1),
('Alexandre Bouchard'   , 9 , 1),
('Catherine Gauthier'   , 10, 1);


INSERT INTO projet.Address (street, apartment, postal_code, province, city) VALUES
('123 Rue Principale'           , 101 , 'H1H 1H1', 'Québec', 'Montréal'),
('456 Avenue Sainte-Catherine'  , 202 , 'H2H 2H2', 'Québec', 'Montréal'),
('789 Rue Sherbrooke'           , 303 , 'H3H 3H3', 'Québec', 'Montréal'),
('101 Rue de la Montagne'       , 404 , 'H4H 4H4', 'Québec', 'Montréal'),
('123 Rue des Fleurs'           , 505 , 'H5H 5H5', 'Québec', 'Montréal'),
('456 Boulevard René-Lévesque'  , 606 , 'H6H 6H6', 'Québec', 'Montréal'),
('789 Avenue du Mont-Royal'     , 707 , 'H7H 7H7', 'Québec', 'Montréal'),
('101 Rue Ontario'              , 808 , 'H8H 8H8', 'Québec', 'Montréal'),
('123 Boulevard Saint-Laurent'  , 909 , 'H9H 9H9', 'Québec', 'Montréal'),
('456 Avenue McGill College'    , 1010, 'H0H 0H0', 'Québec', 'Montréal');


INSERT INTO projet.Borrow (Borrow_DateID, Return_DateID, duration, possible_delay, BookID, SubscriberID, LibraryID) VALUES
(01, 02, 04, 'no', 01, 01, 1),
(02, 03, 10, 'yes', 01, 02, 1),
(03, 04, 11, 'no', 01, 02, 1),
(04, 05, 14, 'no', 03, 02, 1),
(05, 06, 13, 'yes', 05, 05, 1),
(06, 12, 14, 'no', 05, 06, 1),
(07, 08, 06, 'no', 07, 09, 1),
(08, 11, 14, 'no', 09, 09, 1),
(09, 12, 06, 'yes', 10, 09, 1),
(10, 13, 09, 'no', 10, 10, 1);


-- For April (30 days)
INSERT INTO projet.Borrow_Date (day, month, year) VALUES
    (01, 04, 2024),
    (02, 04, 2024),
    (03, 04, 2024),
    (04, 04, 2024),
    (05, 04, 2024),
    (06, 04, 2024),
    (07, 04, 2024),
    (08, 04, 2024),
    (09, 04, 2024),
    (10, 04, 2024),
    (11, 04, 2024),
    (12, 04, 2024),
    (13, 04, 2024),
    (14, 04, 2024),
    (15, 04, 2024),
    (16, 04, 2024),
    (17, 04, 2024),
    (18, 04, 2024),
    (19, 04, 2024),
    (20, 04, 2024),
    (21, 04, 2024),
    (22, 04, 2024),
    (23, 04, 2024),
    (24, 04, 2024),
    (25, 04, 2024),
    (26, 04, 2024),
    (27, 04, 2024),
    (28, 04, 2024),
    (29, 04, 2024),
    (30, 04, 2024),
    (01, 05, 2024),
    (02, 05, 2024),
    (03, 05, 2024),
    (04, 05, 2024),
    (05, 05, 2024),
    (06, 05, 2024),
    (07, 05, 2024),
    (08, 05, 2024),
    (09, 05, 2024),
    (10, 05, 2024),
    (11, 05, 2024),
    (12, 05, 2024),
    (13, 05, 2024),
    (14, 05, 2024),
    (15, 05, 2024),
    (16, 05, 2024),
    (17, 05, 2024),
    (18, 05, 2024),
    (19, 05, 2024),
    (20, 05, 2024),
    (21, 05, 2024),
    (22, 05, 2024),
    (23, 05, 2024),
    (24, 05, 2024),
    (25, 05, 2024),
    (26, 05, 2024),
    (27, 05, 2024),
    (28, 05, 2024),
    (29, 05, 2024),
    (30, 05, 2024),
    (31, 05, 2024),
    (01, 06, 2024),
    (02, 06, 2024),
    (03, 06, 2024),
    (04, 06, 2024),
    (05, 06, 2024),
    (06, 06, 2024),
    (07, 06, 2024),
    (08, 06, 2024),
    (09, 06, 2024),
    (10, 06, 2024),
    (11, 06, 2024),
    (12, 06, 2024),
    (13, 06, 2024),
    (14, 06, 2024),
    (15, 06, 2024),
    (16, 06, 2024),
    (17, 06, 2024),
    (18, 06, 2024),
    (19, 06, 2024),
    (20, 06, 2024),
    (21, 06, 2024),
    (22, 06, 2024),
    (23, 06, 2024),
    (24, 06, 2024),
    (25, 06, 2024),
    (26, 06, 2024),
    (27, 06, 2024),
    (28, 06, 2024),
    (29, 06, 2024),
    (30, 06, 2024);
INSERT INTO projet.Return_Date (day, month, year) VALUES
(01, 04, 2024),
    (02, 04, 2024),
    (03, 04, 2024),
    (04, 04, 2024),
    (05, 04, 2024),
    (06, 04, 2024),
    (07, 04, 2024),
    (08, 04, 2024),
    (09, 04, 2024),
    (10, 04, 2024),
    (11, 04, 2024),
    (12, 04, 2024),
    (13, 04, 2024),
    (14, 04, 2024),
    (15, 04, 2024),
    (16, 04, 2024),
    (17, 04, 2024),
    (18, 04, 2024),
    (19, 04, 2024),
    (20, 04, 2024),
    (21, 04, 2024),
    (22, 04, 2024),
    (23, 04, 2024),
    (24, 04, 2024),
    (25, 04, 2024),
    (26, 04, 2024),
    (27, 04, 2024),
    (28, 04, 2024),
    (29, 04, 2024),
    (30, 04, 2024),
    (01, 05, 2024),
    (02, 05, 2024),
    (03, 05, 2024),
    (04, 05, 2024),
    (05, 05, 2024),
    (06, 05, 2024),
    (07, 05, 2024),
    (08, 05, 2024),
    (09, 05, 2024),
    (10, 05, 2024),
    (11, 05, 2024),
    (12, 05, 2024),
    (13, 05, 2024),
    (14, 05, 2024),
    (15, 05, 2024),
    (16, 05, 2024),
    (17, 05, 2024),
    (18, 05, 2024),
    (19, 05, 2024),
    (20, 05, 2024),
    (21, 05, 2024),
    (22, 05, 2024),
    (23, 05, 2024),
    (24, 05, 2024),
    (25, 05, 2024),
    (26, 05, 2024),
    (27, 05, 2024),
    (28, 05, 2024),
    (29, 05, 2024),
    (30, 05, 2024),
    (31, 05, 2024),
    (01, 06, 2024),
    (02, 06, 2024),
    (03, 06, 2024),
    (04, 06, 2024),
    (05, 06, 2024),
    (06, 06, 2024),
    (07, 06, 2024),
    (08, 06, 2024),
    (09, 06, 2024),
    (10, 06, 2024),
    (11, 06, 2024),
    (12, 06, 2024),
    (13, 06, 2024),
    (14, 06, 2024),
    (15, 06, 2024),
    (16, 06, 2024),
    (17, 06, 2024),
    (18, 06, 2024),
    (19, 06, 2024),
    (20, 06, 2024),
    (21, 06, 2024),
    (22, 06, 2024),
    (23, 06, 2024),
    (24, 06, 2024),
    (25, 06, 2024),
    (26, 06, 2024),
    (27, 06, 2024),
    (28, 06, 2024),
    (29, 06, 2024),
    (30, 06, 2024);

INSERT INTO projet.Order (BookID, SubscriberID, status, LibraryID) VALUES
(1 ,3, 'Fulfilled', 1),
(2 ,3, 'Fulfilled', 1),
(3 ,3, 'Cancelled', 1),
(4 ,3, 'Fulfilled', 1),
(5 ,4, 'Fulfilled', 1),
(6 ,4, 'Fulfilled', 1),
(7 ,2, 'Fulfilled', 1),
(8 ,1, 'Cancelled', 1),
(9 ,6, 'Cancelled', 1),
(10,7, 'Cancelled', 1);



INSERT INTO projet.Library (name) VALUES
('University Library');


commit;
------------------------------------------------------------------------------------
begin;

-- references section
alter table Projet.Subscriber
add constraint fk_subscriber_LibraryID foreign key (LibraryID) references projet.Library(LibraryID),
add constraint fk_subscriber_AddressID foreign key (AddressID) references projet.Address(AddressID);

alter table Projet.Borrow
add constraint fk_borrow_Borrow_DateID foreign key (Borrow_DateID) references projet.Borrow_Date(Borrow_DateID),
add constraint fk_borrow_Return_DateID foreign key (Return_DateID) references projet.Return_Date(Return_DateID),
add constraint fk_borrow_BookID foreign key (BookID) references projet.Books(BookID),
add constraint fk_borrow_SubscriberID foreign key (SubscriberID) references projet.Subscriber(SubscriberID),
add constraint fk_borrow_LibraryID foreign key (LibraryID) references projet.Library(LibraryID);



alter table Projet.Order
add constraint fk_order_SubscriberID foreign key (SubscriberID) references projet.Subscriber(SubscriberID),
add constraint fk_order_LibraryID foreign key (LibraryID) references projet.Library(LibraryID);

commit;
------------------------------------------------------------------------------------
-- add constraints
begin;

ALTER TABLE projet.Borrow
ADD CONSTRAINT check_duration
CHECK (duration <= 14);



ALTER TABLE projet.Borrow
ADD CONSTRAINT check_dates
CHECK (Borrow_DateID < Return_DateID); -- a larger index indicates a larger date value in our case



CREATE OR REPLACE FUNCTION count_fulfilled_orders(subscriber_id INT) RETURNS INT AS $$
DECLARE
    total INT;
BEGIN
    SELECT COUNT(*)
    INTO total
    FROM projet.Order
    WHERE SubscriberID = subscriber_id AND status = 'Fulfilled';

    RETURN total;
END;
$$ LANGUAGE plpgsql;


ALTER TABLE projet.Order
ADD CONSTRAINT max_3_fulfilled_per_subscriber
CHECK (
    status = 'Cancelled' OR
    (
        status = 'Fulfilled' AND
        count_fulfilled_orders(SubscriberID) <= 3
    )
);




commit;

------------------------------------------------------------------------------------